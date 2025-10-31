#!/bin/bash

source /home/caisongbo/anaconda3/bin/activate primer
set -e

# file_list="$1"
root_dir="/s97/lab/woodman/4w/result"
data_dir="/s97/SHARE/splist_nt/nt.p1"
file_dirs="$root_dir/trans_r_nt"
out_dirs="$root_dir/nt_map"

process_file() {
    filename=$(basename "$1")
    #dirname=$(echo "$filename" | awk -F'_' '{print $1"_"$2}')
    #dirname=$(echo "$filename" | cut -c 1-15)
    dirname="${filename%.fa}"
    # echo "Processing $dirname"
    mkdir -p "$out_dirs/$dirname"

    if [ -f $out_dirs/${dirname}/${dirname}_unique.sam ];then
        return
    fi
    
    bowtie2 --mm -p 1 --end-to-end \
    -x $data_dir \
    -f "$1" \
    --no-1mm-upfront \
    --ignore-quals \
    --norc \
    -S $out_dirs/${dirname}/${dirname}.sam \
    --al $out_dirs/${dirname}/${dirname}_aligned.fq \
    --un $out_dirs/${dirname}/${dirname}_unaligned.fq > $out_dirs/${dirname}/${dirname}.log 2>&1

    if [ -f "$out_dirs/${dirname}/${dirname}.sam" ]; then
        :
        # grep "AS:" $out_dirs/${dirname}/${dirname}.sam | grep -v -Ev "(XS:|XM:i:1)" > $out_dirs/${dirname}/${dirname}_unique.sam
        grep "AS:" $out_dirs/${dirname}/${dirname}.sam | grep -v "XS:" > $out_dirs/${dirname}/${dirname}_unique.sam
        rm -f $out_dirs/${dirname}/${dirname}.sam
    else
        echo "Error: SAM file not generated for $dirname" >> error.log
    fi
    
}

export -f process_file
export out_dirs data_dir

# Correctly pass the directory name to process_file
ls "$file_dirs" | parallel --bar -j 50 process_file "$file_dirs/{}"