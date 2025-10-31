import os
import csv
import json
from tqdm.auto import tqdm


def update_csv_data(file_old, file_new, pair_dict):
    with open(file_old, mode='r', encoding='utf-8') as f_in, open(file_new, mode='w', encoding='utf-8', newline='') as f_out:
        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)

        writer.writeheader()
        
        for row in tqdm(reader, desc="Updating csv data"):
            pid = row["X"]
            if pid in pair_dict:
                fields=pair_dict[pid]
                # left_value = fields.get("Left", "F")
                left_value = fields["Left"]
                right_value = fields["Right"]
                if left_value == "T" and right_value == "T":
                    writer.writerow(row)

sub_name = "virus"
work_dir = f"/s97/lab/woodman/{sub_name}/{sub_name}_mate"
file_old = os.path.join(work_dir, f"primer_bank_{sub_name}_v6.6.csv")
file_new = os.path.join(work_dir, f"primer_bank_{sub_name}_v7.csv")

unique_txt="/home/caisongbo/lab_wm/virus/total_virus_id_e7toe3.txt"
# unique_txt="/home/caisongbo/lab_wm/fungi/total_fungi.v2.blast"
# unique_txt="/home/caisongbo/lab_wm/4w/result/total.v3.blast"
results = {}
with open(unique_txt, "r") as f:
    for line in tqdm(f):
        field = line.strip().split("\t")[0]
        infos = field.split("_")
        genome_id = "_".join(infos[:2])
        gene_id = "_".join(infos[2:4])
        seq_type = infos[-4]
        pid="_".join(infos[-3:])  
        #fungi
            # genome_id = "_".join(infos[:2])
            # gene_id = "_".join(infos[2:6])
            # seq_type = infos[6]
            # pid="_".join(infos[7:]) 
        if pid not in results:
            results[pid] = {
                "ID": pid,
                "GENOME_ID": genome_id,
                "GENE_ID": gene_id,
                "Left": "F",
                "Right": "F"
            }

        if seq_type == "LEFT":
            results[pid]["Left"] = "T"
        elif seq_type == "RIGHT":
            results[pid]["Right"] = "T"
print(len(results))

flag=0
for k,v in results.items():
    # if "ORGANISM_NAME" in v:
    if flag<5:
        print(f"Genome_id:{k}")
        print(v)
        flag+=1

update_csv_data(file_old, file_new, results)