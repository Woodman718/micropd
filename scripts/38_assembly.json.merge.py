import csv,os
import glob
import json
from tqdm.auto import tqdm

def txt_to_jsonl(input_txt, output_jsonl, key_field="GENOME_ID"):
    result_dict = {}
    with open(input_txt, 'r') as txt_file, open(output_jsonl, 'w') as jsonl_file:

        for _ in range(2):
            next(txt_file)

        reader = csv.reader(txt_file, delimiter='\t')
        for row in tqdm(reader):
            assembly_accession = row[0]
            taxid = row[5]
            species_taxid = row[6]
            organism_name = row[7]
            strain = row[8]

            if strain == "na":
                strain = "None"
            elif "=" in strain:
                strain = strain.split("=")[1]
            

            json_obj = {
                "GENOME_ID": assembly_accession,
                "TAXID": taxid,
                "SPECIES_TAXID": species_taxid,
                "ORGANISM_NAME": organism_name,
                "STRAIN": strain
            }


            key = json_obj[key_field]
            result_dict[key] = json_obj
            jsonl_file.write(json.dumps(json_obj) + '\n')

    return result_dict

def jsonl_to_dict_by_taxid(jsonl_file_path, key_field="GENOME_ID"):
    data_dict = {}
    
    with open(jsonl_file_path, 'r', encoding='utf-8') as jsonl_file:
        for line in jsonl_file:
            json_obj = json.loads(line)
            taxid = json_obj[key_field]
            data_dict[taxid] = json_obj
    
    return data_dict

def merge_dict_and_write_jsonl(dict1, dict2, output_jsonl):
    merged_dict = dict1.copy()
    merged_dict.update(dict2)

    with open(output_jsonl, 'w') as jsonl_file:
        for key, value in merged_dict.items():
            jsonl_file.write(json.dumps(value) + '\n')

    return merged_dict
  
sub_name = "bacteria"
root_dir = "/s97/lab/woodman"
info_dir = f"{root_dir}/4w/result/info"
input_txt = f"{info_dir}/{sub_name}_assembly_summary_historical.txt"
input_txt2 = f"{info_dir}/{sub_name}_assembly_summary.txt"
output_jsonl = f"{info_dir}/{sub_name}_assembly_GTSS_historical.jsonl"
output_jsonl2 = f"{info_dir}/{sub_name}_assembly_GTSS.jsonl"
output_jsonl3 = f"{info_dir}/{sub_name}_assembly_merge.GTSS.jsonl"

# key_field="TAXID"
result_dict = txt_to_jsonl(input_txt, output_jsonl)
result_dict2 = txt_to_jsonl(input_txt2, output_jsonl2)
print(input_txt,len(result_dict2))
print(input_txt,len(result_dict))
merged_dict = merge_dict_and_write_jsonl(result_dict2, result_dict, output_jsonl3)
len(set(merged_dict.keys()))
