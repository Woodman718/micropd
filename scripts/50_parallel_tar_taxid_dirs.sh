#!/bin/bash

subject="fna"
temp_dirs="/home/caisongbo/lab_wm/4w/result/temp_${subject}_taxid"
gz_dir="/home/caisongbo/lab_wm/4w/result/download/${subject}_taxid"

mkdir -p "$gz_dir"

compress_dir() {
    local temp_dir=$1
    local gz_dir=$2
    local subdir=$(basename "$temp_dir")
    tar -czf "$gz_dir/${subdir}.${subject}.tar.gz" -C "$temp_dir" .
}

export -f compress_dir
export gz_dir subject

find "$temp_dirs" -mindepth 1 -maxdepth 1 -type d | parallel --bar -j 128 compress_dir {} "$gz_dir"

echo "All compression tasks completed."