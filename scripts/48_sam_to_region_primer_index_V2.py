import csv,os
import json
from tqdm.auto import tqdm


csv_file = "./bacteria_mate/primer_bank_bacteria_4w_v3.csv"
tsv_file = "./info/samErrorMatch.tsv"
tsv_file2 = "./info/samNoMatch.tsv"
txt_file = "samErrorMatch.txt"
txt_file2 = "samNoMatch.txt"
output_file = "./bacteria_mate/primer_bank_bacteria_4w_v3.clean.csv"
filter_ids = set()

def read_txt_PDI(txt_file,filter_set):
    with open(txt_file, "r") as f:
        for line in tqdm(f):
            field = line.strip()
            parts = field.split("_")  
            primer_id = "_".join(parts[-3:]) 
            filter_set.add(primer_id)
    return filter_set

def read_tsv_PDI(tsv_file,filter_set):
    with open(tsv_file, "r") as tsv:
        reader = csv.reader(tsv, delimiter="\t")
        next(reader) 
        for row in tqdm(reader):
            if row:
                primer_id = '_'.join(row[0].split('_')[-3:])
                filter_set.add(primer_id)
    return filter_set
    
read_txt_PDI(txt_file,filter_ids)
read_txt_PDI(txt_file2,filter_ids)

with open(csv_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    header = next(reader)
    writer.writerow(header)
    
    for row in tqdm(reader):
        if row[0] not in filter_ids:
            writer.writerow(row)
            
#####
root_dir = "/s97/lab/woodman/4w/result"
original_csv = f'{root_dir}/bacteria_mate/primer_bank_bacteria_4w_v3.clean.csv'
output_dir = f'{root_dir}/primer/csv2fasta'
total_fasta = f"{root_dir}/total.v3.fasta"
total_jsonl = f"{root_dir}/total.v3.jsonl"

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

genome_data = {}
with open(total_fasta, 'w') as fasta_file, open(total_jsonl, 'w') as jsonl_file:
    # Read and process the CSV file
    with open(original_csv, 'r') as original_file:
        reader = csv.DictReader(original_file)
        # Organize data by GENOME_ID
        for row in tqdm(reader, desc="Processing CSV rows"):
            genome_id = row['GENOME_ID']
            gene_id = row['GENE_ID']
            left_seq = row['PRIMER_LEFT_SEQUENCE']
            right_seq = row['PRIMER_RIGHT_SEQUENCE']
            seq_id = row['X']

            # FASTA
            left_fasta = f">{genome_id}_{gene_id}_LEFT_{seq_id}\n{left_seq}"
            right_fasta = f">{genome_id}_{gene_id}_RIGHT_{seq_id}\n{right_seq}"
            
            fasta_file.write(left_fasta + '\n')
            fasta_file.write(right_fasta + '\n')
            
            # Create a list of entries for each GENOME_ID
            if genome_id not in genome_data:
                genome_data[genome_id] = []
            genome_data[genome_id].extend([left_fasta, right_fasta])
            
            json_entry = {
                "PID": seq_id,
                "LEFT_SEQ": left_seq,
                "RIGHT_SEQ": right_seq,
                "GENOME_ID": genome_id,
                "GENE_ID": gene_id
            }
            jsonl_file.write(json.dumps(json_entry) + '\n')
            
# Write each GENOME_ID group to a separate FASTA file
for genome_id, sequences in tqdm(genome_data.items(), desc="Write FASTA"):
    fasta_filename = f"{output_dir}/{genome_id}.fa"
    with open(fasta_filename, 'w') as individual_fasta_file:
        individual_fasta_file.write('\n'.join(sequences) + '\n')
        