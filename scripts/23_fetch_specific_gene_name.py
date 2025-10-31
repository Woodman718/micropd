import os
import pandas as pd
import argparse

ap = argparse.ArgumentParser(description="Process Fastq files in a specified directory.")
ap.add_argument('-f', '--txt', type=str, required=True, help='The txt directory.')
ap.add_argument('-s', '--sam', type=str, default=None, help='The sam directory.')
ap.add_argument('-o', '--output', type=str, default=None, help='The output directory.')
args = vars(ap.parse_args())

input_txt = args["txt"]
sam_dir = args["sam"]
output_dir = args["output"]

genome_out = f"{output_dir}/specific_genome_bacteria.txt"
txt_f = os.path.basename(input_txt)
genome_name = txt_f.split('.txt')[0]
genome_id = genome_name[:15]

# Check if input txt file is empty
if os.path.getsize(input_txt) == 0:
    with open(genome_out, 'a') as f:
        f.write(f"{genome_name} is empty\n")
    exit(1)

specific_gene_id = f"{output_dir}/{genome_id}_specific_gene_id.txt"
unique_sam_path = os.path.join(sam_dir, genome_id, f"{genome_id}_unique.sam")

# Check if unique.sam file exists and is not empty
if os.path.getsize(unique_sam_path) > 0:
    # Read gene data from the input txt file
    gen_read = pd.read_csv(input_txt, sep="\t", header=None, index_col=None)
    # Read unique.sam data
    uni_read = pd.read_table(unique_sam_path, header=None, sep="\t", skip_blank_lines=True).reset_index(drop=True)
    uni_read_name = uni_read[0].str.split("/", expand=True)[0]
    uni_read_name_table = uni_read_name.value_counts().reset_index()
    uni_read_name_table.columns = ['uni_read_name', 'Freq']

    # Identify specific genes
    specific_gene = set(gen_read[0]).intersection(set(uni_read_name_table['uni_read_name']))
    specific_gene_read = gen_read[gen_read[0].isin(specific_gene)].cop()
    specific_gene_read.columns = ['geneid', 'Tread']

    # Map unique read names to frequencies and add Vread column
    specific_gene_vread = uni_read_name_table[uni_read_name_table['uni_read_name'].isin(specific_gene)]
    specific_gene_vread = specific_gene_vread.set_index('uni_read_name')['Freq']
    specific_gene_read['Vread'] = specific_gene_read['geneid'].map(specific_gene_vread)

    # Filter rows with matching Tread and Vread counts
    specific_gene_read = specific_gene_read.dropna()
    specific_gene_read['T_Vread'] = specific_gene_read['Tread'] - specific_gene_read['Vread']
    specific_gene_read_efficient = specific_gene_read[specific_gene_read['T_Vread'] == 0]

    # If no specific gene data, log as empty and exit with code 1
    if specific_gene_read_efficient.empty:
        with open(genome_out, 'a') as f:
            f.write(f"{genome_name} specific_gene_read_efficient is empty\n")
        exit(1)
    else:
        # Write specific gene data to output file
        specific_gene_name = specific_gene_read_efficient['geneid']
        specific_gene_name.to_csv(specific_gene_id, sep="\t", index=False, header=False)
else:
    with open(genome_out, 'a') as f:
        f.write(f"{genome_name}_unique.sam is empty\n")
    exit(1)
