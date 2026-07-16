import subprocess
import sys
import os

configfile: "config.yaml"

# =========================
# ENV VARIABLES
# =========================
BED = config["bed_file"]
RESULTS = config["results"]
SAMPLES_LIST = config["samples_list"]

def get_seq_list(wildcards):
    """
    This function is called AFTER the checkpoint completes.
    It reads the sequence list created by reverse_fasta.
    
    This must be defined BEFORE it's used in rule all,
    but it executes AFTER the checkpoint finishes.
    """
    checkpoint_output = checkpoints.reverse_fasta.get().output.seq_list
    
    with open(checkpoint_output, "r") as f:
        seq_ids = [line.strip() for line in f if line.strip()]
    
    return seq_ids
    
rule all:
    input: f"{RESULTS}/data/astral/species_tree.newick"
    #input: expand(f"{RESULTS}/data/trees/{{seq_id}}.treefile", seq_id=get_seq_list)
    
# ---- rule to generate index files ----
rule make_index:
    input:
        vcf = config["vcf"],
        genome = config["genome"]
    output:
        vcf_csi = config["vcf"] + ".csi",
        vcf_tbi = config["vcf"] + ".tbi",
        genome_dict = config["genome"].replace(".fasta", "") + ".dict"
    log:
        f"{RESULTS}/logs/make_index.log"
    conda:
        "envs/vcf_tools.yaml"

    resources:
        mem_mb = 100000,
        cpus = 16,
        runtime = 300
    shell:
        """
        bcftools index {input.vcf} 
        gatk IndexFeatureFile -I {input.vcf} 
        gatk CreateSequenceDictionary -R {input.genome} 2> {log}
        """

if config["filter"]:

	rule gatk_filter_vcf:
		input:
			vcf = config["vcf"],
			genome = config["genome"]
		output:
			vcf_filtered = f"{RESULTS}/data/VCF_DP.vcf.gz",
		conda:
			"envs/vcf_tools.yaml"
		params:
			DP_min = config["DP_min"],
			DP_max = config["DP_max"]
		log:
			f"{RESULTS}/logs/VCF_DP.log"
		resources:
			mem_mb = 300000,
			cpus = 16,
			runtime = 200
		shell:
			"""
			gatk --java-options "-Xmx10g" VariantFiltration \
			-R {input.genome} \
			-V {input.vcf} \
			-O {output.vcf_filtered} \
			--genotype-filter-name "DPFILTER" \
			--genotype-filter-expression "DP<{params.DP_min} ||  DP>{params.DP_max}" \
			--set-filtered-genotype-to-no-call true \
			--verbosity INFO 2> {log} 
			"""
			
	rule index_DPfilter:
		input:
			vcf_filtered = f"{RESULTS}/data/VCF_DP.vcf.gz"
		output:
			vcf_index_bcf = f"{RESULTS}/data/VCF_DP.vcf.gz.csi",
		conda:
			"envs/vcf_tools.yaml"
		resources:
			mem_mb = 300000,
			cpus = 1,
			runtime = 60
		shell:
			"""
			bcftools index {input.vcf_filtered}
			"""
			
	rule filter_NA_vcf:
		input:
			vcf_filtered = f"{RESULTS}/data/VCF_DP.vcf.gz"
		output:
			vcf_filteringBase = f"{RESULTS}/data/VCF_DP_NA.vcf.gz",
			vcf_index_csi = f"{RESULTS}/data/VCF_DP_NA.vcf.gz.csi"
		params:
			Na_rate = config["Na_rate"]
		conda:
			"envs/vcf_tools.yaml"
		log:
			f"{RESULTS}/logs/VCF_Filters_DP_NA.log"
		resources:
			mem_mb = 10000,
			cpus = 4,
			runtime = 200
		shell:
			"""
			bcftools view \
			-i 'F_MISSING<{params.Na_rate}' \
			-O z \
			{input.vcf_filtered} \
			-o {output.vcf_filteringBase} 2> {log}
			
			bcftools index {output.vcf_filteringBase}
			"""
			
	FINAL_VCF = rules.filter_NA_vcf.output.vcf_filteringBase

else:
    FINAL_VCF = config["vcf"]

rule get_consensus:
    wildcard_constraints:
        sample = "|".join(SAMPLES_LIST)  # Force same samples
    input:
        vcf = FINAL_VCF,
        bed = config["bed_file"],
        genome = config["genome"],
        fai = config["fai"]
    output:
        fasta = f"{RESULTS}/data/by_samples/{{sample}}.fasta"
    params:
        haplotype = config["haplotype"]
    conda:
        "envs/vcf_tools.yaml"
    log:
        f"{RESULTS}/logs/by_sample/{{sample}}.log"
    resources:
        mem_mb = 10000,
        cpus = 2,
        runtime = 120
    shell:
        """
        set -euo pipefail
        samtools faidx {input.genome} -r {input.bed} | 
        bcftools consensus \
        --absent "N"\
        --haplotype {params.haplotype} \
        --sample {wildcards.sample} \
        -M "N" \
        {input.vcf} > {output.fasta} 2> {log}
        """

rule compute_Nper:
    input:
        fasta = f"{RESULTS}/data/by_samples/{{sample}}.fasta"
    output: 
        statout = f"{RESULTS}/data/by_samples/{{sample}}_Nper.txt"
    log:
        f"{RESULTS}/logs/by_sample/{{sample}}_Nper.log"
    conda:
        "envs/phylo_tools.yaml"
    resources:
        mem_mb = 5000,
        cpus = 1,
        runtime = 180
    shell:
        """
        python -c "from fasta_tools.utils import get_N_percent; get_N_percent('{input.fasta}', '{output.statout}', '{wildcards.sample}')"
        """

rule merge_Nper :
    input:
        Nper_table = expand(f"{RESULTS}/data/by_samples/{{sample}}_Nper.txt", sample=SAMPLES_LIST)
    output:
        Nper_tablout = f"{RESULTS}/data/by_samples/ALL_Nper.txt"
    conda:
        "envs/phylo_tools.yaml"
    log:
        f"{RESULTS}/logs/ALL_Nper.log"
    resources:
        mem_mb = 5000,
        cpus = 1,
        runtime = 30
    shell:
        """
        python << 'EOF'
import pandas as pd

files = "{input.Nper_table}".split()
df_list = [pd.read_csv(file, sep='\t') for file in files]

df_merged = df_list[0]
for df in df_list[1:]:
    df_merged = df_merged.merge(df, on='seq_id', how='outer')

df_merged.to_csv('{output.Nper_tablout}', sep='\t', index=False)
EOF
        """

rule Filter_Nper:
    wildcard_constraints:
        sample = "|".join(SAMPLES_LIST)  # Force same samples
    input:
        Nper_table = f"{RESULTS}/data/by_samples/ALL_Nper.txt",
        fasta_list = f"{RESULTS}/data/by_samples/{{sample}}.fasta"
    output:
        fasta_filtered = f"{RESULTS}/data/by_samples/filtered_{{sample}}.fasta"
    conda:
        "envs/phylo_tools.yaml"
    params:
        Nper_threshold = config["max_na"],
        min_taxa = config["min_taxa"]
    log:
        f"{RESULTS}/logs/{{sample}}_fasta_filtered.log"
    resources:
        mem_mb = 5000,
        cpus = 1,
        runtime = 30
    shell:
        """
        python << 'EOF'
import pandas as pd
from fasta_tools.utils import Select_Seq

tablin = pd.read_csv('{input.Nper_table}',  sep="\t", index_col=0)
tablin_tag=tablin.where(tablin < {params.Nper_threshold}, other=pd.NA)
sum_ind = tablin_tag.transpose().count()
perc_filtered_seq = sum_ind*100 / len(tablin.columns)
Seq_list=perc_filtered_seq[perc_filtered_seq > {params.min_taxa}].index.tolist()
Select_Seq('{input.fasta_list}', Seq_list, '{output.fasta_filtered}')
EOF
        """

### Checkpoint to reorganize fasta by sequence ID   ################

checkpoint reverse_fasta:
    wildcard_constraints:
        sample = "|".join(SAMPLES_LIST)
    input:
        fasta_files = expand(f"{RESULTS}/data/by_samples/filtered_{{sample}}.fasta", sample=SAMPLES_LIST)
    output:
        done = f"{RESULTS}/data/by_seq/reverse_fasta.done",
        seq_list = f"{RESULTS}/data/by_seq/seq_id_list.txt"
    log:
        f"{RESULTS}/logs/reorganize_by_seq.log"

    resources:
        mem_mb = 10000,
        cpus = 1,
        runtime = 160
    conda:
        "envs/phylo_tools.yaml"
    shell:
        """
        python script/reverse.py \
            {RESULTS}/data/by_seq \
            {input.fasta_files}

        touch {output.done}
        """
 
rule multiple_alignment:
    input:
        fasta = f"{RESULTS}/data/by_seq/filtered_rev_{{seq_id}}.fasta"
    output:
        aligned_fasta = f"{RESULTS}/data/align/{{seq_id}}.fasta"
    log:
        f"{RESULTS}/logs/alignment/{{seq_id}}.log"
    conda:
        "envs/phylo_tools.yaml"
    threads: 16
    resources:
        mem_mb = 50000,
        cpus = 16,
        runtime = 240
    shell:
        """
        clustalo -i {input.fasta} -o {output.aligned_fasta} -t DNA --outfmt=fasta --threads={threads} --force 2> {log}
        """

rule tree_by_seq:
    input:
        aligned_fasta = f"{RESULTS}/data/align/{{seq_id}}.fasta"
    output:
        tree = f"{RESULTS}/data/trees/{{seq_id}}.treefile"
    params:
        prefix = f"{RESULTS}/data/trees/{{seq_id}}"
    log:
        f"{RESULTS}/logs/trees/{{seq_id}}.log"
    conda:
        "envs/phylo_tools.yaml"
    threads: 32
    resources:
        mem_mb = 50000,
        cpus = 32,
        runtime = 300
    shell:
        """

        iqtree -s {input.aligned_fasta} --redo-tree -m MFP -nt {threads} -pre {params.prefix} 2> {log}
        """

rule astral_tree:
    input:
        trees = expand(f"{RESULTS}/data/trees/{{seq_id}}.treefile", seq_id=get_seq_list)
    output:
        astral_tree = f"{RESULTS}/data/astral/species_tree.tre",
        newick_tree = f"{RESULTS}/data/astral/species_tree.newick"
    params:
        mappingfile = config["mapping_file"]
    conda:
        "envs/phylo_tools.yaml"
    log:
        f"{RESULTS}/logs/astral_tree.log"
    resources:
        mem_mb = 50000,
        cpus = 4,
        runtime = 240
    shell:
        """
        cat {input.trees} > {output.newick_tree}
        astral -i {output.newick_tree} -a {params.mappingfile} -o {output.astral_tree} 2> {log}
        """
