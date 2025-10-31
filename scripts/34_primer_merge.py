import os
import glob
import pandas as pd
from tqdm import tqdm


root_dir="/home/caisongbo/s97/lab/woodman/4w/result"
data_dir=f"{root_dir}/primer/ann_uniref"
file_name="primer_bank_bacteria_4w_v1.csv"
# data_dir=f"{root_dir}/demo/test"
# file_name="4w_t3.csv"

output_file=os.path.join(root_dir,file_name)


file_list = glob.glob(f"{data_dir}/*.csv")
# file_list = glob.iglob(f"{data_dir}/*.csv")
# file_list = os.listdir(data_dir)

id_prefix = "PRIMER_PAIR_"
primer_bank = []
tmp=0
batch_size = 1e6 
#total=len(file_list)
for f_path in tqdm(file_list, desc="Merge files"):
    primer_info = pd.read_csv(f_path)
    primer_info['KINGDOM'] = primer_info['KINGDOM'].fillna('NA')

    num_records = len(primer_info)
    primer_info.insert(0, 'X', [f"{id_prefix}{i}" for i in range(tmp, tmp + num_records)])
    tmp += num_records
    primer_bank.append(primer_info)

    if len(primer_bank) >= batch_size:
        pd.concat(primer_bank, ignore_index=True).to_csv(
            output_file, mode='a', header=not os.path.exists(output_file), index=False
        )
        primer_bank = []

if primer_bank:
    pd.concat(primer_bank, ignore_index=True).to_csv(
        output_file, mode='a', header=not os.path.exists(output_file), index=False
    )