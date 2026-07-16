## to do
 - fill the profile file to add all slurm parameters by rule (allow the pipeline to run outside slurm if needed)
 - add containers 

# Pipeline_VCF2Tree

A Snakemake 8 pipeline that turns a multi-sample VCF into a species phylogeny, designed to run on the IFB cluster (SLURM).

Starting from a VCF + reference genome (and .fai index) + a BED file of target regions, the pipeline builds per-sample, per-region consensus sequences, filters them, aligns them, builds a gene tree per region, and summarizes all gene trees into a single species tree with ASTRAL.

---

## 1. What the Snakefile does

The pipeline runs the following steps, in order:

1. **`make_index`** — Indexes the input VCF (`bcftools index`, `gatk IndexFeatureFile`) and creates the genome sequence dictionary (`gatk CreateSequenceDictionary`).

2. **`gatk_filter_vcf`** *(optional, only if `filter: true` in `config.yaml`)* — Filters genotypes by read depth using `gatk VariantFiltration`. Positions with `DP < DP_min` or `DP > DP_max` are set to no-call. Output: `VCF_DP.vcf.gz`.

3. **`index_DPfilter`** — Indexes the depth-filtered VCF.

4. **`filter_NA_vcf`** — Removes variant sites with too much missing data across samples, using `bcftools view -i 'F_MISSING<Na_rate'`. Output: `VCF_DP_NA.vcf.gz`.
   *(If `filter: false`, steps 2–4 are skipped and the original VCF is used directly.)*

5. **`get_consensus`** — For each sample, extracts the target regions from the genome (`samtools faidx` + `bed_file`) and builds a consensus FASTA sequence with `bcftools consensus`, using the chosen `haplotype` mode. Missing genotypes become `N`.

6. **`compute_Nper`** — Computes the percentage of `N` (missing) bases per sequence, per sample.

7. **`merge_Nper`** — Merges all per-sample `%N` tables into a single table (`ALL_Nper.txt`).

8. **`Filter_Nper`** — Filters out sequences that are too incomplete: a sequence must have less than `max_na`% missing data, and be present in at least `min_taxa`% of samples to be kept.

9. **`reverse_fasta`** (checkpoint) — Reorganizes the filtered, per-sample FASTA files into per-sequence-ID FASTA files (one file per genomic region, containing all samples) via `script/reverse.py`. This is a Snakemake *checkpoint* because the list of surviving sequence IDs is only known after filtering.

10. **`multiple_alignment`** — Aligns each per-region FASTA across samples with `clustalo`.

11. **`tree_by_seq`** — Builds a gene tree for each aligned region with `iqtree` (`-m MFP`, automatic model selection).

12. **`astral_tree`** — Concatenates all gene trees and infers the final species tree with `astral`, using the `mapping_file` (sample → species/population mapping).

**Final output:**
```
results/data/astral/species_tree.newick
```

### Conda environments
Two environments are used, defined per rule:
- `envs/vcf_tools.yaml` — `bcftools`, `samtools`, `gatk4`
- `envs/phylo_tools.yaml` — `clustalo`, `iqtree`, `astral`, `pandas`, and the custom `fasta_tools` module

---

## 2. How to run the test

The repository ships with a small test dataset in `input/data_test/` (test VCF, test genome, BUSCO-derived BED file, mapping file) and `config.yaml` is already pointed at these files by default, so you can run the pipeline as-is with no configuration changes.

### On the IFB cluster (SLURM)
```bash
sbatch run_skm.sh
```
`run_skm.sh` loads the `snakemake` and `python` modules, sets up a local conda package cache, and launches Snakemake with:
```bash
CONDARC=.condarc snakemake --use-conda --conda-frontend mamba --executor slurm --profile profiles/slurm
```
It uses the SLURM executor plugin together with the cluster profile in `profiles/slurm/`, so each rule is submitted as its own SLURM job using the `resources:` (mem, cpus, runtime) declared in the Snakefile.

### Locally / on a single machine (not available for now because all rules come with slurm parameters directly in the snakefile)
```bash
snakemake --use-conda --conda-frontend mamba --cores <N>
```

### Useful checks before running for real
```bash
# Dry run: shows which rules/jobs would execute
snakemake -n

# Visualize the DAG
snakemake --dag | dot -Tpng > dag.png
```

Results of the test run land in `results/` (created automatically, or see the pre-computed example in `run_test/results/`). Logs for each rule are written to `results/logs/`.

---

## 3. How to run your own data

### Step 1 — Add your files to the `input/` folder
Create a subfolder (e.g. `input/my_data/`) and place:

| File | Description |
|---|---|
| VCF (`.vcf.gz`) | Multi-sample variant file, bgzipped |
| Genome (`.fasta`) | Reference genome the VCF was called against |
| Genome index (`.fasta.fai`) | `samtools faidx genome.fasta` |
| BED file | 4 columns: `chr  start  end  name` — the regions/loci you want a gene tree for |
| Mapping file | Sample → species/population mapping, used by ASTRAL |

### Step 2 — Edit `config.yaml`

```yaml
# paths
vcf: "input/my_data/my_variants.vcf.gz"
genome: "input/my_data/my_genome.fasta"
fai: "input/my_data/my_genome.fasta.fai"
bed_file: "input/my_data/my_regions.bed"
results: "results"
samples_list: ["sample1", "sample2", "sample3", ...]   # must match sample names in the VCF
mapping_file: "input/my_data/my_mapping_file"

# VCF filtering
filter: false          # set to true to enable the DP-based genotype filter
DP_min: 10              # min depth to keep a genotype (below → N)
DP_max: 150             # max depth to keep a genotype (above → N)
BaseQual: 30             # min base quality
Na_rate: 0.5             # max fraction of missing genotypes allowed per site

# Consensus sequence parameters
haplotype: "I"     # 1 = first allele, 2 = second allele, R = REF, A = ALT, I = IUPAC code

# Sequence filtering (post-consensus)
max_na: 20          # % missing max allowed per sequence
min_taxa: 95         # % of samples that must share a sequence to keep it

# Hyper parameters
minperc_sample_by_seq: 95.0
```

Key points:
- **`samples_list`** must exactly match the sample names in your VCF header — this drives the per-sample wildcards throughout the pipeline.
- **`filter`** toggles the optional depth-based genotype filtering step (rules 2–4 above). Leave it `false` if your VCF is already filtered.
- **`haplotype`** controls how heterozygous/ambiguous genotypes are encoded in the consensus FASTA (see comments in `config.yaml` for the 5 available modes).
- **`max_na` / `min_taxa`** control how aggressively incomplete sequences are dropped before alignment — raise `min_taxa` for a stricter, more complete matrix; lower it to keep more loci at the cost of more missing data.

### Step 3 — Run
Same commands as the test run (Section 2), either via `sbatch run_skm.sh` on IFB, or directly with `snakemake --use-conda --cores <N>` locally.

### Step 4 — Check the output
The final species tree is written to:
```
results/data/astral/species_tree.newick
```
Intermediate outputs (per-sample consensus FASTAs, alignments, gene trees, logs) are kept under `results/data/` and `results/logs/` for inspection/debugging.
