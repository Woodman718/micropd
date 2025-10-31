import pysam
import os
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file from fasta")
ap.add_argument("-r", "--background", required=True, help="input txt file include gene names")
ap.add_argument("-o", "--output_dir", required=True, help="directory to store output fasta files")
ap.add_argument("-L", "--log_dir", required=True, help="directory to store log files")

args = vars(ap.parse_args())

input_fa = args["input"]
input_txt = args["background"]
output_dir = args["output_dir"]
log_idr = args["log_dir"]

genome_id= os.path.basename(input_fa)
log_file_path = os.path.join(log_idr,f"{genome_id}.txt")
with open(input_txt, 'r') as fi_txt:
    gene_ids = [line.strip() for line in fi_txt]

with pysam.FastxFile(input_fa, "r") as fi_fa:
    for entry in fi_fa:
        sequence_name = entry.name
        sequence = entry.sequence

        if sequence_name in gene_ids:
            output_file_path = os.path.join(output_dir, f"{sequence_name}.fasta")
            with open(output_file_path, "w") as fo:
                fo.write(f">{sequence_name}\n")
                fo.write(f"{sequence}\n")

log_message = f"{genome_id} Fasta files have been saved.\n"
with open(log_file_path, "w") as log_file:
    log_file.write(log_message)
    
print(log_message)