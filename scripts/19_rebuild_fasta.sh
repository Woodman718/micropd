#!/bin/bash

# Load conda environment
source /home/caisongbo/s1/anaconda3/bin/activate web

set -e 
# Define paths
root_dir="/home/caisongbo/lexicmap730"
merged_fa_dir="$root_dir/bacteria_cds/output/cds_hit_filter"  # Directory with merged fasta files
# merged_fa_dir="$root_dir/bacteria_cds/test"
gene_name_dir="$root_dir/unique/gene_name"  # Directory with gene name txt files
output_dir="$root_dir/remove_redundancy_fa"  # Directory for output genome-specific fasta files
py_script="$root_dir/merge2fasta.py"  # Python script path
 

# Run the Python script
python "$py_script" -m "$merged_fa_dir" -g "$gene_name_dir" -o "$output_dir"