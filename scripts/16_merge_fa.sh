#!/bin/bash

# Define directory paths
root_dir="/home/caisongbo/lexicmap730"
input_dir="$root_dir/bacteria_cds/4w_cds_filter_fa"
output_dir="$root_dir/4w_cds_merge_fa"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Calculate the total number of files
total_files=$(ls $input_dir | wc -l)
files_per_batch=1000
batch_number=1

# Initialize file count
file_count=0

# Create a temporary file to store file names
temp_file_list="file_list.txt"

# Clear the temporary file (if it exists)
> $temp_file_list

# Iterate over all files in the input directory
for file in $input_dir/*; do
    echo $file >> $temp_file_list
    ((file_count++))

    # Merge every 1000 files
    if [ $file_count -eq $files_per_batch ]; then
        cat $(cat $temp_file_list) > $output_dir/merged_batch_${batch_number}.fa
        ((batch_number++))
        file_count=0
        > $temp_file_list
    fi
done

# Process any remaining files
if [ $file_count -gt 0 ]; then
    cat $(cat $temp_file_list) > $output_dir/merged_batch_${batch_number}.fa
fi

# Remove the temporary file
rm $temp_file_list

echo "File merging complete, total $batch_number batches."