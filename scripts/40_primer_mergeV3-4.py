import csv,os
import json
from tqdm.auto import tqdm
from collections import Counter
import matplotlib.pyplot as plt

def read_txt_PDI(txt_file):
    filter_set = set()
    with open(txt_file, "r") as f:
        for line in tqdm(f):
            field = line.strip()
            parts = field.split("_") 
            primer_id = "_".join(parts[-3:]) 
            filter_set.add(primer_id)
    return filter_set
  
root_dir = "/home/caisongbo/s97/lab/woodman/virus"
txt_file = f"{root_dir}/virus.v3.e7.blast"
unique_txt = f"{root_dir}/total_virus_id_e7toe3.txt"
thresholds = 3 #default:1
csv_file = f"{root_dir}/virus_mate/primer_bank_virus_v3.csv"
output_file = f"{root_dir}/virus_mate/primer_bank_virus_v4.csv"

with open(txt_file, 'r') as f:
    lines = f.readlines()

fields = []
for line in lines:
    parts = line.strip().split("\t")
    fields.append(parts[0])

line_counts = Counter(fields)
unique_lines = [line for line, count in line_counts.items() if count <= thresholds]

with open(unique_txt, 'w') as file:
    for line in unique_lines:
        file.write(f"{line}\n")
        
print(f"Unique lines based on first two columns: {len(unique_lines)}({len(lines)})")

unique_ids = read_txt_PDI(unique_txt )
with open(csv_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    header = next(reader)
    writer.writerow(header)
    
    for row in tqdm(reader):
        # if row[0] not in filter_ids: 
        if row[0] in unique_ids:
            writer.writerow(row)