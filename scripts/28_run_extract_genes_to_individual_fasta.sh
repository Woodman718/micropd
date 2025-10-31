#!/bin/bash


root_dir="/s3/woodman/Primer"
fasta_dir="$root_dir/speci_genome_fa"
txt_dir="$root_dir/speci_genome_id"
out_fa_dir="$root_dir/speci_gene_fa"
log_dir="$root_dir/logs"
py_file="/home/caisongbo/s3/Primer/extract_genes_to_individual_fasta.py"

process_file() {
    filename=$(basename "$1")
    genemo_id=$(echo "$filename" | cut -c 1-15)
    txt_file="$txt_dir/${genemo_id}_specific_gene_id.txt"
    log_file="${log_dir}/${genemo_id}.fasta.txt"

    # Check whether "txt_file" exists and the corresponding fasta file does not exist
    if [ ! -f "$log_file" ]; then
        python "$py_file" -i "$1" \
                        -r "$txt_file" \
                        -o "$out_fa_dir" \
                        -L "$log_dir"
    else
        echo "$genemo_id has been processed."
    fi
}

mkdir -p "$out_fa_dir"

export -f process_file
export root_dir fasta_dir txt_dir out_fa_dir log_dir py_file 

ls "$fasta_dir" | parallel --bar -j 64 process_file "$fasta_dir/{}"