#!/bin/bash

set -e

# work documentation
root_dir="/s97/lab/woodman/4w/result"

csv_name="$root_dir/all_csv_nameUniref.txt"
csv_dir="$root_dir/primer/ann_uniref"
output_dir="$root_dir/primer/score"
int_med_dir="$root_dir/primer/intermediate_score"
py_file="/home/caisongbo/opt/score_primers.py"


process_file() {
    csv_path="$1"
    file_name=$(basename "${csv_path%_primer3.csv}")
    output_file="$output_dir/${file_name}_score.csv"
    if [ -f "$output_file" ]; then
        echo "$(date '+%Y-%m-%d') - ${file_name} - existed"
    else
        python "$py_file" "$csv_path" "$output_dir" "$int_med_dir" "$file_name"
    fi  
}

mkdir -p $output_dir
mkdir -p $int_med_dir

export -f process_file
export root_dir csv_dir output_dir py_file int_med_dir


cat "$csv_name" | xargs -P 218 -I {} bash -c 'process_file "$csv_dir/{}"'