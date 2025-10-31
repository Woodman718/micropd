import os
import pandas as pd
from tqdm import tqdm


root_dir="/home/caisongbo/lab_wm/4w/result"
# data_dir=f"{root_dir}/demo/csv"
data_dir=f"{root_dir}/primer/bacteria"

# output_dir=f"{root_dir}/demo"
file_name="bacteria_bank.csv"
output_file=os.path.join(root_dir,file_name)

file_list = os.listdir(data_dir)
primer_bank = []
ids,tmp=0,0
flag = True

for f in tqdm(file_list, desc="Merge files"):
    f_path=os.path.join(data_dir,f)
    if os.path.getsize(f_path) > 3:
        info_dat = f.split("_") 
        genomeid = f"{info_dat[0]}_{info_dat[1]}"
        geneid = f"{info_dat[2]}_{info_dat[3]}" 
        
        dat = pd.read_csv(f_path, index_col=0)
        
        primer_info = dat.T.reset_index(drop=True) 
        
        if not flag:
            ids += len(primer_info)
            primer_info.index = [f"PRIMER_PAIR_{i}" for i in range(tmp,ids)]
            tmp+=len(primer_info)
        primer_info['GENOME_ID'] = genomeid 
        primer_info['GENE_ID'] = geneid 
        
        
        primer_bank.append(primer_info)

primer_bank = pd.concat(primer_bank, ignore_index=flag)

primer_bank.to_csv(output_file) 