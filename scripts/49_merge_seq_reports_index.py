import csv
import os
import json
import pickle
from tqdm.auto import tqdm
from collections import defaultdict


def generate_and_merge_reports(reports_dir, filter_dir, work_dir):
    merged_data = {}
    temp_data = {}
    
    file_list = os.listdir(reports_dir)
    for filename in tqdm(file_list,total=len(file_list)):
        Genome_id = filename[:15]
        temp_path = os.path.join(filter_dir, f"{Genome_id}.json")

        if filename.endswith("_seq_report.jsonl"):
            seq_report_file = os.path.join(reports_dir, filename)
            with open(seq_report_file, 'r') as f:
                for line in f:
                    json_obj = json.loads(line)
                    genbank_accession = json_obj['genbankAccession']
                    assembly_accession = json_obj['assemblyAccession']
                    refseq_accession = json_obj['refseqAccession']
                    MoleculeType = json_obj['assignedMoleculeLocationType']
                    # 检查是否已经存在key
                    if genbank_accession in merged_data:
                        print(f"Warning: Duplicate genbankAccession '{genbank_accession}' found in file {filename}. Existing entry will be overwritten.")
                    merged_data[genbank_accession] = [assembly_accession, refseq_accession,MoleculeType]
                    temp_data[genbank_accession] = [assembly_accession, refseq_accession,MoleculeType]
            with open(temp_path, 'w') as tempf:
                json.dump(temp_data, tempf, indent=4)
                
    merged_report_file = os.path.join(work_dir, "bacteria_seq_report.json")
    with open(merged_report_file, 'w') as jsonf:
        json.dump(merged_data, jsonf, indent=4)
    
    report_file = os.path.join(work_dir, "bacteria_seq_report.pk")
    with open(report_file, "wb") as pkf:
        pickle.dump(dict(merged_data), pkf)


if __name__ == "__main__":
    root_dir="/s97/lab/woodman/4w/result"
    reports_dir = os.path.join(root_dir, "sequence_reports")
    filter_dir = os.path.join(root_dir, "filter_seq_report")
    os.makedirs(filter_dir, exist_ok=True)
    
    generate_and_merge_reports(reports_dir, filter_dir, root_dir)
    print("Complete!")