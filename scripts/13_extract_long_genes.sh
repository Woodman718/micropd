#!/bin/bash

# Define directories and output paths
root_dir="/s1/store91/woodman/jupyter/tzd/Prokka"
fasta_dir="${root_dir}/bacteria_cds/4w"
output="$root_dir/gene_len/4w"

# Export variables for parallel processing
export root_dir fasta_dir output

# Print the start time
echo "$(date +%F_%H:%M:%S)"

# Create the output directory if it doesn't exist
mkdir -p "$output"

# Parallel processing of files
ls "$fasta_dir" | parallel -j --bar 32 '
    i="{}"
    gene_name="${i%.fasta}"

    # Check if the output file already exists
    if [ ! -f "$output/$gene_name.txt" ]; then
        # Convert FASTA to tabular format with seqkit
        seqkit fx2tab -n -i -l "$fasta_dir/$i" > "$output/$gene_name.txt"
    else
        echo "$gene_name exists."
    fi
'

# Print the end time
echo "$(date +%F_%H:%M:%S)"
echo "Merging files..."

# Merge the output files into a single file
find "$output" -type f -name "*.txt" | xargs -n 1000 cat >> ../4w_cds_gene_lengths.txt

# Extract gene IDs of sequences longer than 10,000 bases
awk '$2 > 10000 {print $1}' ../4w_cds_gene_lengths.txt > ../4w_cds_gene_more_than_1w.txt

echo "Finished."