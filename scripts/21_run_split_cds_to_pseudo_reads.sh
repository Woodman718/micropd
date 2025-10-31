#!/bin/bash

# Activate the specified Anaconda environment `primer`
source /home/caisongbo/anaconda3/bin/activate primer

set -e 

# Define directories and variables
root_dir="/home/caisongbo"
swd_dir="$root_dir/s97/SHARE/caisongbo"

# Path to the Python script for processing files
bins="$root_dir/opt/split_cds_to_pseudo_reads.py"
data_dir="$swd_dir/remove_redundancy_fa/4w"
fastq_dir="$swd_dir/result/split150/bacteria/fastq"
txt_dir="$swd_dir/result/split150/bacteria/txt"

# Get the current date and time for logging
target_date="$(date +%F)"
current="$(date +%H%M%S)"
log_dir="$swd_dir/result/logs/${target_date}"
log_file="$log_dir/${current}_process.log"
exist_file="$log_dir/${current}_exist.log"
works="40"

# Create necessary directories
mkdir -p $fastq_dir
mkdir -p $txt_dir
mkdir -p $log_dir

# Export variables for use in the parallel environment
export data_dir fastq_dir txt_dir bins exist_file log_file

# List all files in the data directory and process each file in parallel using GNU Parallel
ls "$data_dir"| parallel  --bar -j ${works} '
    j="{}"
    i=${j%.fasta}   # Remove `.fasta` suffix to get the base filename
    result_file="$fastq_dir/${i}.fastq"  # Path to the target fastq file

    # If the result file does not exist, perform the conversion
    if [ ! -e "$result_file" ]; then
        echo $i >> $log_file 

        # Run the Python script for file conversion
        python ${bins} -i $data_dir/$j -o $result_file -j $txt_dir/${i}.txt

    else
        echo $result_file >> $exist_file
    fi
'

# After the script completes, write a finish message to the log file
echo "Script execution finished." >> $log_file