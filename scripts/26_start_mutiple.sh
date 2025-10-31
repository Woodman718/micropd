#!/bin/bash
target_date=$(date +%F)
log_file="${target_date}_start.log"
exist_file="${target_date}_exist.log"
start_count=0
file_list="$1"
bin_file="batch_diamond_uniref.sh"
result_dir="./result/4w"
batch_size=50
batch=()

mkdir -p $result_dir

> "$log_file"
> "$exist_file"

for j in $(cat ${file_list}); do
    file_dir="$result_dir/$j"
    if [ ! -d "$file_dir" ]; then
        batch+=($j)
        if [ ${#batch[@]} -eq $batch_size ]; then
	    tmp_list_file=$(mktemp -p ./diamond)
            printf "%s\n" "${batch[@]}" > $tmp_list_file
            sbatch $bin_file $tmp_list_file
            echo "${batch[@]}" >> "$log_file"
            start_count=$((start_count + batch_size))
            batch=()
        fi
    else
        echo $file_dir >> $exist_file
    fi
done

# If there are still unsubmitted entries remaining
if [ ${#batch[@]} -gt 0 ]; then
    tmp_list_file=$(mktemp -p ./diamond)
    printf "%s\n" "${batch[@]}" > $tmp_list_file
    sbatch $bin_file $tmp_list_file

    echo "${batch[@]}" >> "$log_file"
    start_count=$((start_count + ${#batch[@]}))
fi

sed -i "1i\Total directories started: $start_count"  "$log_file"
sed -i "1i\Operation completed at $(date +%F_%H:%M:%S)"  "$log_file"

head -n 2 "$log_file"