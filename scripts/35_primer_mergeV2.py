import os
import glob
import pandas as pd
from tqdm import tqdm

root_dir = "/home/caisongbo/s97/lab/woodman/4w/result"
gc_len_dir = f"{root_dir}/speci_gene_gc_len"
data_dir=f"{root_dir}/primer/score"
file_name="primer_bank_bacteria_4w_v2.csv"
output_file = os.path.join(root_dir, file_name)
id_prefix = "PRIMER_PAIR_"
primer_bank = []
tmp = 0
batch_size = 1e5 

#file_list = glob.glob(f"{data_dir}/*.csv")
file_list = os.listdir(data_dir)
#file_list.sort()

gene_gc_len_files = glob.glob(f"{gc_len_dir}/*_specific_gene_gc_len.txt")
gene_gc_len_dict = {}

for gc_len_file in gene_gc_len_files:
    df_gc_len = pd.read_csv(
        gc_len_file, sep=r"\s+", header=None, names=["gene_id", "GENE_LENGTH", "GENE_GC"]
    )
    gene_gc_len_dict.update(df_gc_len.set_index("gene_id").to_dict(orient="index"))


for f_path in tqdm(file_list, desc="Merge score files"):
    primer_info = pd.read_csv(f_path)
    primer_info['KINGDOM'] = primer_info['KINGDOM'].fillna('NA')

    if 'score' in primer_info.columns:
        primer_info.rename(columns={'score': 'SCORE'}, inplace=True)
    
    # Obtain the length and GC content of the gene based on the GENE_ID
    primer_info["GENE_LENGTH"] = primer_info["GENE_ID"].map(
        lambda gene_id: gene_gc_len_dict.get(gene_id, {}).get("GENE_LENGTH", "NA")
    )
    primer_info["GENE_GC"] = primer_info["GENE_ID"].map(
        lambda gene_id: gene_gc_len_dict.get(gene_id, {}).get("GENE_GC", "NA")
    )

    num_records = len(primer_info)
    primer_info.insert(0, 'X', [f"{id_prefix}{i}" for i in range(tmp, tmp + num_records)])
    tmp += num_records
    primer_bank.append(primer_info)

    if len(primer_bank) >= batch_size:
        pd.concat(primer_bank, ignore_index=True).to_csv(
            output_file, mode='a', header=not os.path.exists(output_file), index=False
        )
        primer_bank = []

# Write the remaining data
if primer_bank:
    pd.concat(primer_bank, ignore_index=True).to_csv(
        output_file, mode='a', header=not os.path.exists(output_file), index=False
    )