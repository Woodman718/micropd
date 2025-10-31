#!/bin/bash
#SBATCH --job-name=uniref_multiple
#SBATCH --output=./info/diamond_%j.out
#SBATCH --error=./info/diamond_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=7
#SBATCH --mem=64G

##version：1.4
source /public/home/caisongbo/anaconda3/bin/activate prokka_diamond
set -e

list_file="$1"
job_id=$SLURM_JOB_ID
bin_dir="./uniref_annotator"
in_dir="./faa_dir/4"
out_dir="./result/4w"
target_date="$(date +%F)"
current="$(date +%H%M%S)"
log_dir="./logs/${target_date}"
exist_file="$log_dir/${job_id}_exist.log"
log_file="$log_dir/${job_id}.log"

db_dir="./dmnd_db"
uniref90="$db_dir/uniref90.dmnd"
uniref50="$db_dir/uniref50.dmnd"
works="7"

export bin_dir in_dir out_dir
export db_dir uniref90 uniref50
export works exist_file log_file


# Check if diamond is installed and executable
if ! command -v diamond &> /dev/null; then
    echo "diamond could not be found"
    exit 1
fi

if [[ ! -x "$(command -v diamond)" ]]; then
    echo "diamond is not executable"
    exit 1
fi

mkdir -p $log_dir
#mkdir -p $out_dir
# Get job details and log them

parallel -j 2 '
i="{}"
result_dir="$out_dir/$i"

if [ ! -d "$result_dir" ]; then
    mkdir -p $result_dir
    echo $i >> $log_file

  python ${bin_dir}/uniref_annotator.py \
         ${in_dir}/$i/$i.faa \
         --seqtype prot \
         --uniref90db ${uniref90} \
         --uniref50db ${uniref50} \
         --diamond-options "--threads ${works}" \
         --temp $result_dir \
         --out $result_dir/$i.annotated
else
    echo $result_dir >> $exist_file
fi
' < "$list_file"

scontrol show job $job_id > "$log_file"
echo "Script execution finished." >> $log_file
