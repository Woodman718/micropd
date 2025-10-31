#!/bin/bash

set -e

# Activate the conda environment
source /home/caisongbo/anaconda3/bin/activate cdhit

# Define variables
#target_date="$(date +%F)"
target_date="2024-09-28"
root_dir="/s98/lab/woodman/primer/cdhit"
input_dir="$root_dir/12w_cds_merge_fa"
output_dir="$root_dir/output_12w/$target_date"

# Create base output and filter directories
mkdir -p "$output_dir"

# Set parameters for cd-hit-est and seqkit
identity=0.95            # -c parameter
memory=400000            # -M parameter (in MB)
threads=24               # -T parameter
word_length=10           # -n parameter
alignment=0.9            # -aS parameter

# Loop through each file in the input directory
for fast_file in "$input_dir"/*; do
    echo `date +%H%M%S`
    gene_name=$(basename "${fast_file%.fa}")
    out_file="$output_dir/${gene_name}.fasta"
    logs="$output_dir/${gene_name}.log"

    if [[ -f "$out_file" ]]; then
        echo "Output file '$out_file' already exists. Skipping."
        continue  # Skip to the next iteration of the loop
    fi
    # Run cd-hit-est
    cd-hit-est -i "${fast_file}" \
               -o "${out_file}" \
               -G 0 \
               -n ${word_length} \
               -aS ${alignment} \
               -c ${identity} \
               -M ${memory} \
               -d 0 \
               -r 1 \
               -g 1 \
               -T ${threads} \
               -B 0 > "$logs" 2>&1
              # &> "$logs"
done