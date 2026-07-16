import os
import sys
import fasta_tools

# -----------------------
# Arguments
# -----------------------

output_dir = sys.argv[1]
fasta_files = sys.argv[2:]

# -----------------------
# Create output dir
# -----------------------

os.makedirs(output_dir, exist_ok=True)

# -----------------------
# Store seq by ID
# -----------------------

seq_dict = {}

for fasta_file in fasta_files:

    basename = os.path.basename(fasta_file)

    current_sample = (
        basename
        .replace('filtered_', '')
        .replace('.fasta', '')
    )

    print(f"Processing sample: {current_sample}")

    dico_fasta = fasta_tools.fasta_dict(fasta_file)

    for seq_id, sequence in dico_fasta.items():

        if seq_id not in seq_dict:
            seq_dict[seq_id] = []

        seq_dict[seq_id].append((current_sample, sequence))

print(f"Total unique sequences: {len(seq_dict)}")

# -----------------------
# Write fasta by sequence
# -----------------------

for seq_id, samples_data in seq_dict.items():

    seq_file=f"{output_dir}/seq_id_list.txt"
    with open(seq_file, "a") as seq_out:
        seq_out.write(seq_id + "\n")

    output_file = f"{output_dir}/filtered_rev_{seq_id}.fasta"
    with open(output_file, "w") as out:

        for sample_name, sequence in samples_data:

            header = f"{sample_name}"

            out.write(">"+header + "\n")
            out.write(sequence + "\n")