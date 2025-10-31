#!/bin/bash

source /home/caisongbo/anaconda3/bin/activate primer
set -e

files="$1"
root_dir="/home/caisongbo"
data_dir="$root_dir/s98/SHARE/bowtie2_index/bacteria"
file_dirs="$root_dir/lab_wm/4w/result/split150/bacteria/fastq"
out_dirs="$root_dir/lab_wm/4w/result/mapping/Bacteria"

process_file() {
    filename=$(basename "$1")
    #dirname=$(echo "$filename" | awk -F'_' '{print $1"_"$2}')
    genemo_id=$(echo "$filename" | cut -c 1-15)
    # echo "Processing $dirname"
    mkdir -p "$out_dirs/$dirname"
		
    if [ -f $out_dirs/${dirname}/${dirname}_unique.sam ];then
        return
    fi
    
    # -f , process fasta files
    # --sam-append-comment
    # Append FASTA/FASTQ comment to SAM record, where a comment is everything after the first space in the read name.
    bowtie2 --mm -p 1 --end-to-end \
    -x $data_dir \
    -U "$1" \
    -S $out_dirs/${dirname}/${dirname}.sam \
    --al $out_dirs/${dirname}/${dirname}_aligned.fq \
    --un $out_dirs/${dirname}/${dirname}_unaligned.fq > $out_dirs/${dirname}/${dirname}.log 2>&1

    if [ -f "$out_dirs/${dirname}/${dirname}.sam" ]; then
        grep "AS:" $out_dirs/${dirname}/${dirname}.sam | grep -v "XS:" > $out_dirs/${dirname}/${dirname}_unique.sam
    else
        echo "Error: SAM file not generated for $dirname" >> error.log
    fi
    rm -f $out_dirs/${dirname}/${dirname}.sam
}

export -f process_file
export out_dirs data_dir

# Correctly pass the directory name to process_file
ls "$file_dirs" | parallel --bar -j 100 process_file "$file_dirs/{}"