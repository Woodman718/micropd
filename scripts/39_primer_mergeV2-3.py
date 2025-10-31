import csv,os
import json
from tqdm.auto import tqdm

def load_gene_mapping(jsonl_file):
    gene_mapping = {}
    with open(jsonl_file, 'r') as infile:
        for line in infile:
            record = json.loads(line)
            gene_id = record["GENE_ID"]
            gene_name = record.get("GENE_NAME", "")
            gene_location = record.get("GENE_LOCATION", "")
            gene_mapping[gene_id] = (gene_name, gene_location)
    return gene_mapping
    
 def update_primer_csv(main_table, jsonl_file, output_file):

    gene_mapping = load_gene_mapping(jsonl_file)

    with open(main_table, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        header = next(reader)
        gene_id_index = header.index("GENE_ID")
        new_header = (
            header[:gene_id_index + 1]
            + ["GENE_NAME", "GENE_LOCATION"]
            + header[gene_id_index + 1:]
        )
        writer.writerow(new_header)

        for row in tqdm(reader, desc="Updating CSV file"):
            gene_id = row[gene_id_index]
            gene_name, gene_location = gene_mapping.get(gene_id, ("", ""))
            new_row = (
                row[:gene_id_index + 1]
                + [gene_name, gene_location]
                + row[gene_id_index + 1:]
            )
            writer.writerow(new_row)
   
work_dir="/s97/lab/woodman/virus/virus_mate"
main_table = os.path.join(work_dir,"primer_bank_virus_v2.csv")
jsonl_file = f"{work_dir}/virus_gff.jsonl"
output_file = os.path.join(work_dir,"primer_bank_virus_v3.csv")

update_primer_csv(main_table, jsonl_file, output_file)