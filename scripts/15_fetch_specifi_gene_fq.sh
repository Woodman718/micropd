#!/bin/bash

# Define directories and output paths
root_dir="/home/caisongbo/lexicmap730"
out_fa_dir="$root_dir/bacteria_cds/4w_cds_filter_fa"
fil_fa_dir="$root_dir/bacteria_cds/4_cds_1w_plus"
py_file="$root_dir/script/fetch_specifi_gene_fq.py"
fasta_dir="/s3/SHARE/woodman/Prokka2/bacteria_cds/4w"
genome_txt="$root_dir/fasta_39901"

# Create output directories if they don't exist
mkdir -p "$out_fa_dir"
mkdir -p "$fil_fa_dir"

# Export variables for parallel processing
export root_dir out_fa_dir fil_fa_dir py_file fasta_dir

# Print the start time
echo "$(date +%F_%H:%M:%S)"

# Parallel processing of files
cat "$genome_txt" | parallel --bar -j 32 '
    i="{}"
    gene_name="${i%.fasta}"

    # Check if the output file already exists
    if [ ! -f "$out_fa_dir/$gene_name.fa" ]; then
        # Run the Python script to process the gene file
        python "$py_file" -i "$fasta_dir/${gene_name}.fasta" \
                           -r "$root_dir/4w_cds_gene_more_than_1w.txt" \
                           -o "$out_fa_dir/${gene_name}.fa" \
                           -f "$fil_fa_dir/${gene_name}.fa"
    else
        echo "$gene_name already processed."
    fi
'

# Print the end time
echo "$(date +%F_%H:%M:%S)"
echo "Finished."