#!/bin/bash

source /home/caisongbo/anaconda3/bin/activate primer3
set -e 

target_date="$(date +%F)"
# root_dir="/home/caisongbo/lab_wm/4w/result"
# 12w
root_dir="/s97/SHARE/caisongbo/result"
fasta_dir="$root_dir/speci_genome_fa"
out_dir="$root_dir/primer/bacteria"
log_dir="$root_dir/logs/$target_date"
py_file="/home/caisongbo/opt/primer_design_primer3.py"


process_file() {
    # filename=$(basename "$1")
    # genemo_id="${filename%.fasta}"
    # genemo_id=$(echo "$filename" | cut -c 1-15)
    # out_fasta="${out_fa_dir}/${genemo_id}.fasta"
    # fil_fasta="${fil_fa_dir}/${genemo_id}.fasta"
    python $py_file -i "$1" -o ${out_dir} --log $log_dir

}

mkdir -p "$out_dir"
mkdir -p "$log_dir"

export -f process_file
export root_dir fasta_dir out_dir log_dir py_file 


ls "$fasta_dir" | parallel --bar -j 128 process_file "$fasta_dir/{}"