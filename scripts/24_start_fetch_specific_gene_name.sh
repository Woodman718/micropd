#!/bin/bash

set -e

root_dir="/home/caisongbo/lab_wm/4w/result"
txt_dir="$root_dir/split150/bacteria/txt"
data_dir="$root_dir/mapping/Bacteria"
out_dir="$root_dir/speci_genome_id"

py_file="/home/caisongbo/scripts/fetch_specific_gene_name.py"

process_file() {
    filename=$(basename "$1")
    gene_name="${filename%.txt}"
    genemo_id=$(echo "$filename" | cut -c 1-15)
    sam_dir="$data_dir/$genemo_id"

    if [ -f $out_dir/${genemo_id}_specific_gene_id.txt ]; then
        return
    fi

    # Run Python script and check exit status
    if python "$py_file" -f "$1" -s "$data_dir" -o "$out_dir"; then
        # If Python script is successful, log gene_name
        echo "$gene_name" >> "${root_dir}/genome_specific_id_bacteria.txt"
    fi
}

mkdir -p $out_dir

export -f process_file
export data_dir root_dir out_dir py_file

ls "$txt_dir" | parallel --bar -j 64 process_file "$txt_dir/{}"