import csv,os
import glob
import json
from tqdm.auto import tqdm

def jsonl_to_dict_by_taxid(jsonl_file_path, key_field="GENOME_ID"):
    data_dict = {}
    
    with open(jsonl_file_path, 'r', encoding='utf-8') as jsonl_file:
        for line in jsonl_file:
            json_obj = json.loads(line)
            taxid = json_obj[key_field]
            data_dict[taxid] = json_obj
    
    return data_dict

def update_primer(main_table, assembly_dict, output_file):
    with open(main_table, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # 读取表头并调整位置
        header = next(reader)
        genome_id_index = header.index("GENOME_ID")
        taxid_id_index = header.index("TAXID")
        loc_index = taxid_id_index + 1
        new_header = (
            header[:loc_index]
            + ["SPECIES_TAXID", "STRAIN"]
            + header[loc_index:]
        )
        writer.writerow(new_header)

        # 逐行处理
        for row in tqdm(reader, desc="Updating CSV file"):
            genome_id = row[genome_id_index]
            if genome_id in assembly_dict:
                genome_data = assembly_dict[genome_id]
                species_taxid = genome_data.get("SPECIES_TAXID", "None")
                strain = genome_data.get("STRAIN", "None")
            new_row = (
                row[:loc_index]
                + [species_taxid, strain]
                + row[loc_index:]
            )
            writer.writerow(new_row)

sub_name = "bacteria"
root_dir = "/s97/lab/woodman"
info_dir = f"{root_dir}/4w/result/info"
input_txt = f"{info_dir}/{sub_name}_assembly_summary_historical.txt"
input_txt2 = f"{info_dir}/{sub_name}_assembly_summary.txt"
output_jsonl = f"{info_dir}/{sub_name}_assembly_GTSS_historical.jsonl"
output_jsonl2 = f"{info_dir}/{sub_name}_assembly_GTSS.jsonl"
output_jsonl3 = f"{info_dir}/{sub_name}_assembly_merge.GTSS.jsonl"

work_dir = f"{root_dir}/4w/result"
file_v5 = f"{work_dir}/{sub_name}_mate/primer_bank_{sub_name}_4w_v5.csv"
updated_file=f"{work_dir}/{sub_name}_mate/primer_bank_{sub_name}_4w_v6.csv"

assembly_dict = jsonl_to_dict_by_taxid(output_jsonl3)
update_primer(file_v5, assembly_dict, updated_file)