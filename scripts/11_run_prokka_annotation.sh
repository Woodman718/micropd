#!/bin/bash
# conda activate Bio prokka_diamond
output='./output'
list_file="/s2/woodman/Prokka2/labels/part_00_aj"

cd $output
echo $(date +%F_%H:%M:%S)
parallel -j 8 '
    name=$(basename "{}" | sed 's/_genomic.*//')
    #mkdir -p "$name"
    prokka "/s2/woodman/Prokka2/twelve_bacterial_fa_qc/{}" --prefix "$name" --cpus 8 --kingdom Bacteria
' < "$list_file"

echo $(date +%F_%H:%M:%S)