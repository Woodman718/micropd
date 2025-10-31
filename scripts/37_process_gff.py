import os
import glob
import json
from tqdm.auto import tqdm

def parse_gff_to_jsonl(gff_dir, jsonl_file):
    with open(jsonl_file, 'w') as outfile:
        gff_files = glob.glob(os.path.join(gff_dir, "*.gff"))
        for gff_file in tqdm(gff_files, desc="Parsing GFF files", total=len(gff_files)):
            genome_id = os.path.basename(gff_file)[:15]
            with open(gff_file, 'r') as infile:
                for line in infile:
                    if line.startswith(">"):
                        break
                    if line.startswith("#"):
                        continue
                    
                    columns = line.strip().split('\t')
                    chrom = columns[0]
                    start = columns[3]
                    end = columns[4]
                    attributes = columns[8]
                    
                    gene_id = None
                    gene_name = None 
                    for attribute in attributes.split(';'):
                        if '=' in attribute: 
                            key, value = attribute.split('=')
                            if key == "ID":
                                gene_id = value
                            elif key == "Name":
                                gene_name = value
                    
                    if gene_id:
                        if gene_name is None:
                            gene_name = gene_id
                        gene_location = f"{chrom}:{start}..{end}"
                        record = {
                            "GENOME_ID": genome_id,
                            "GENE_ID": gene_id,
                            "GENE_NAME": gene_name,
                            "GENE_LOCATION": gene_location
                        }
                        outfile.write(json.dumps(record) + '\n')
                        
work_dir="/s97/lab/woodman/virus/virus_mate"
gff_dir = os.path.join(work_dir,"genome_browser")
jsonl_file = f"{work_dir}/virus_gff.jsonl"
parse_gff_to_jsonl(gff_dir, jsonl_file)