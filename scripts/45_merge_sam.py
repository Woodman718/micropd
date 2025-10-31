import os
import csv
import pandas as pd
from tqdm.auto import tqdm

root_dir="/s97/lab/woodman/4w/result/demo"
dir_a = f"{root_dir}/nt_A"
dir_b = f"{root_dir}/nt_B"
output_dir = f"{root_dir}/merged_sam"
os.makedirs(output_dir, exist_ok=True)


sam_files_a = sorted(os.listdir(dir_a))
sam_files_b = sorted(os.listdir(dir_b))

def read_sam(filepath):
    rows = []
    with open(filepath, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows

for sam_file_a, sam_file_b in tqdm(zip(sam_files_a, sam_files_b)):
    if sam_file_a[:15] == sam_file_b[:15]:
        base_name = sam_file_a[:15]
        output_file = os.path.join(output_dir, f"{base_name}.sam")

        rows_a = read_sam(os.path.join(dir_a, sam_file_a))
        rows_b = read_sam(os.path.join(dir_b, sam_file_b))

        df_a = pd.DataFrame(rows_a)
        df_b = pd.DataFrame(rows_b)
        merged_df = pd.concat([df_a, df_b], ignore_index=True)

        merged_df = merged_df.drop_duplicates(subset=0, keep=False)
        merged_df = merged_df[~merged_df.apply(lambda x: "XS:" in "\t".join(map(str, x)), axis=1)]
        # merged_df = merged_df[~merged_df.apply(lambda x: "XS:" in x.to_string(), axis=1)]
        # merged_df = merged_df.sort_values(by=0, ascending=True)
        merged_df.sort_values(by=0, ascending=True, inplace=True)
        merged_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"The merge is complete and the result is saved in the directory: {output_dir}.")