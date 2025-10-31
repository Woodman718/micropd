import csv
import os,glob
import json,pickle
import pandas as pd
import numpy as np
import gzip,zipfile,tarfile
import shutil
from tqdm.auto import tqdm
import pysam


def validate_table(input_csv, output_csv, accession_taxid_map, output_file):
    genome_id_dict = {}
    with open(input_csv, "r") as csv_in, open(output_csv, "w", newline="") as csv_out:
        reader = csv.DictReader(csv_in)
        writer = csv.DictWriter(csv_out, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in tqdm(reader):
            genome_id = row['GENOME_ID']
            primer_id = row['X']
            taxid=row[subject]
            if not taxid:
                try : 
                    if genome_id in accession_taxid_map:
                        row[subject] = accession_taxid_map[genome_id]
                        writer.writerow(row)
                        if genome_id not in genome_id_dict:
                            genome_id_dict[genome_id] = set()
                        genome_id_dict[genome_id].add(primer_id)
                    else:
                        print(f"The[{primer_id}]: {genome_id} not in accession_taxid_map.")
                except Exception as e:
                    print(f"Error processing {primer_id}: {genome_id}, error: {e}")
            else:
                writer.writerow(row)
                
    if genome_id_dict:
        with open(output_file, "w") as json_file:
            json.dump({k: list(v) for k, v in genome_id_dict.items()}, json_file, indent=4)

    return genome_id_dict
    
def process_table(file_path):
    taxid_dict = {}
    # data = pd.read_csv(file_path)
    # for _, row in data.iterrows():
    with open(file_path, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in tqdm(reader):
            taxid = row['TAXID']
            genome_id = row['GENOME_ID']
            if taxid not in taxid_dict:
                taxid_dict[taxid] = set()
            taxid_dict[taxid].add(genome_id)

    filtered_taxids = {taxid: genome_ids for taxid, genome_ids in taxid_dict.items() if len(genome_ids) >= 2}

    if filtered_taxids:
        with open(outfile, "w") as json_file:
            json.dump({k: list(v) for k, v in filtered_taxids.items()}, json_file, indent=4)

    return filtered_taxids

subject = "TAXID"
input_csv = "/home/caisongbo/lab_wm/4w/result/bacteria_mate/primer_bank_bacteria_4w_v4.5.csv"
output_csv= "/home/caisongbo/lab_wm/4w/result/bacteria_mate/primer_bank_bacteria_4w_v5.5.csv"
output_file = f"/home/caisongbo/lab_wm/4w/result/genome_{subject}_empty.json"
cache_path = "/home/caisongbo/lab_wm/4w/result/info/bac_assembly_lineage_historical.pkl"

with open(cache_path, 'rb') as f:
    cache_data = pickle.load(f)
    
accession_taxid_map = cache_data['accession_taxid_map']
lineage_dat = cache_data['lineage_dat']

#filtered_taxids = process_table(file_path)
filtered_taxids = validate_table(input_csv, output_csv, accession_taxid_map, output_file)
print("TGENOME_IDs:", filtered_taxids)

