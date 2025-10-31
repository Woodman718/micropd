import csv
import os,glob
import json,pickle
import pandas as pd
import numpy as np
import shutil
from tqdm.auto import tqdm
import pysam

def make_outdir(dir_path):
    os.makedirs(dir_path,exist_ok=True)
    
def split_table_by_organism(input_file, output_dir):
		make_outdir(output_dir)
    with open(input_file, mode='r') as infile:
        reader = csv.DictReader(infile)

        file_handles = {}

        try:
            for row in tqdm(reader):
                organism_name = row["ORGANISM_NAME"]
                taxid = row["TAXID"]
                #if taxid != "2740746":
                #    continue
                #else:
                # 如果没有对应的文件句柄，创建新文件
                if taxid not in file_handles:
                    output_file = os.path.join(output_dir, f"{taxid}.csv")
                    file_handle = open(output_file, mode='w', newline='')
                    writer = csv.DictWriter(file_handle, fieldnames=reader.fieldnames)
                    writer.writeheader()
                    file_handles[taxid] = (file_handle, writer)
                
                # 写入当前行到对应的文件
                file_handles[taxid][1].writerow(row)

        finally:
            # 关闭所有打开的文件句柄
            for file_handle, _ in tqdm(file_handles.values()):
                file_handle.close()

def process_csv(file_path,target_columns,flag="None"):
    taxid_dict = {}
    genome_dict = {}
    gene_dict = {}

    with open(file_path, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in tqdm(reader):
            taxid = row[target_columns[0]]
            if flag == "fungi":
                gene_id = row[target_columns[3]]
            else:
                gene_id = row[target_columns[4]]
            genome_id = row[target_columns[2]]

            if taxid not in taxid_dict:
                taxid_dict[taxid] = set()
            if genome_id not in genome_dict:
                genome_dict[genome_id] = set()
            if gene_id not in gene_dict:
                gene_dict[gene_id] = set()

            taxid_dict[taxid].add(genome_id)
            genome_dict[genome_id].add(gene_id)
            gene_dict[gene_id].add(genome_id)

    return taxid_dict, genome_dict,gene_dict

def split_fasta(input_fa,gene_to_taxid,root_dir,sub_name="faa"):
    # 存储 TAXID 对应的输出句柄
    taxid_file_handles = {}
    work_dir = f"{root_dir}/4w/result/temp_{sub_name}_taxid"
    temp_dir = f"{root_dir}/4w/result/{sub_name}_taxid"
    make_outdir(work_dir)
    print(work_dir)
    # 使用 pysam.FastxFile 解析输入 fasta 文件
    with pysam.FastxFile(input_fa, "r") as fi_fa:
        for entry in tqdm(fi_fa):
            sequence_name = entry.name
            sequence = entry.sequence

            # 检查是否为需要提取的基因
            if sequence_name in gene_to_taxid:
                taxid = gene_to_taxid[sequence_name]
                if taxid not in taxid_file_handles:
                    # 打开一个新的 TAXID 文件
                    if sub_name == "faa":
                        taxid_file = f"{taxid}.faa"
                    else:
                        taxid_file = f"{taxid}.fasta"
                    taxid_file_path = os.path.join(work_dir, taxid_file)
                    taxid_file_handles[taxid] = open(taxid_file_path, "w")
                
                # 写入对应的 TAXID 文件
                taxid_file_handles[taxid].write(f">{sequence_name}\n{sequence}\n")

    # 关闭所有打开的文件
    for handle in taxid_file_handles.values():
        handle.close()

    print("基因序列提取完成！")
    return work_dir,temp_dir

def parse_gtf_line(line, file_prefix):
    fields = line.strip().split("\t")

    if "gene" not in fields:
        return None
    gene_index = fields.index("gene")

    Region = fields[0].split(" ")[0]
    start = fields[gene_index + 1]
    end = fields[gene_index + 2]
    icon = fields[gene_index + 4] 
    gene_id = fields[gene_index + 6].split("gene_id")[1].split('"')[1]  # 第7列，去掉引号和分号
    loc = f"{file_prefix}_{gene_id}"
    record = f"{Region}\t{start}\t{end}\t{loc}\t{loc}\t{icon}\n"

    return record

def parse_gff_line(line):
    if line.startswith(">"):
        return "stop"
    if line.startswith("#"):
        return "skip"
    columns = line.strip().split('\t')
    chrom = columns[0]
    start = columns[3]
    end = columns[4]
    ico = columns[6]
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
    
    # if gene_id:
    if gene_name is None:
        gene_name = gene_id
    
    record = f"{chrom}\t{start}\t{end}\t{gene_id}\t{gene_name}\t{ico}\n"

    return record

def parse_gff_to_txt(gff_dir,out_dir,json_data,suffix="gff"):
    gff_files = os.listdir(gff_dir)
    file_handles = {}
    make_outdir(out_dir)
    for genome_id in tqdm(gff_files, desc=f"Parsing {suffix} files", total=len(gff_files)):
        if genome_id not in json_data:
            continue

        gff_file = os.path.join(gff_dir,genome_id,f"{genome_id}.{suffix}")
        bed_path = os.path.join(out_dir, f"{genome_id}.bed")

        with open(gff_file, 'r') as infile, open(bed_path, 'w') as outfile:
            for line in infile:
                if suffix == "gff":
                    record = parse_gff_line(line)
                elif suffix == "gtf":
                    record = parse_gtf_line(line, genome_id)

                if record:
                    if record == "stop":
                        break
                    elif record == "skip":
                        continue
                    else:
                        outfile.write(record)

def mv_file(json_data,temp_dir,work_dir,missing_dict,suffix=".fasta",flag='None'):
    base_dir = os.path.basename(work_dir)
    for taxid, gff_names in tqdm(json_data.items(),total=len(json_data)):
        # 为每个 taxid 创建临时目录
        taxid_temp_dir = os.path.join(temp_dir, taxid)
        os.makedirs(taxid_temp_dir, exist_ok=True)
        
        for gff_name in gff_names:
            if base_dir == flag:
                gff_file_path = os.path.join(work_dir, gff_name,f'{gff_name}{suffix}')
            else:
                gff_file_path = os.path.join(work_dir, f'{gff_name}{suffix}')
                
            dest_file_path = os.path.join(taxid_temp_dir, f'{gff_name}{suffix}')
            try:
                # 复制 GFF 文件到临时目录
                shutil.copy(gff_file_path, dest_file_path)
            except FileNotFoundError:
                # 记录缺失的 GFF 文件
                missing_dict.setdefault(taxid, set()).add(gff_name)

    if missing_dict:
        print("Missing files:")
        for k, vs in missing_dict.items():
            print(f'TaxID {k}: {", ".join(vs)}')
    else:
        print("All files processed successfully.")

subject = "bacteria"
target_columns  = ["TAXID","ORGANISM_NAME","GENOME_ID","CDS_ID", "GENE_ID","GENE_GC"]
root_dir = "/s97/lab/woodman"
csv_dir=f"{root_dir}/4w/result/download/csv_taxid"

genome_browser_dir = f"{root_dir}/4w/result/{subject}_mate/genome_browser"
temp_fna_dir = f"{root_dir}/4w/result/temp_fna_taxid"
temp_gff_dir = f"{root_dir}/4w/result/temp_gff_taxid"
temp_bed_dir = f"{root_dir}/4w/result/temp_bed_taxid"
bed_dir = f"{root_dir}/4w/result/bed_taxid"
val_csv=f"{root_dir}/4w/result/bacteria_mate/primer_bank_bacteria_4w_v6.csv"
input_pro_fa=f"{root_dir}/4w/result/{subject}_mate/4w_specificGene_all_{subject}_pro.fa"
input_fa=f"{root_dir}/4w/result/{subject}_mate/4w_specificGene_all_{subject}.fa"
print(val_csv,f"\n",input_fa)

missing_gff = {}
missing_faa = {}
missing_fna = {}
missing_fa = {}
missing_bed = {}

taxid_dict, genome_dict, gene_dict = process_csv(val_csv,target_columns)
gene_dict2 = {k: list(v)[0] for k, v in gene_dict.items()}
json_data = {k: list(v) for k, v in taxid_dict.items()}

#csv
split_table_by_organism(val_csv, csv_dir)
#fna
mv_file(json_data,temp_fna_dir,genome_browser_dir,missing_fna,suffix=".fna",flag='genome_browser')

#gtf or gff
mv_file(json_data,temp_gff_dir,genome_browser_dir,missing_gff,suffix=".gtf",flag='genome_browser')

#faa
work_dir,temp_dir = split_fasta(input_pro_fa,gene_dict2,root_dir,sub_name="faa")
mv_file(json_data,temp_dir,work_dir,missing_faa,suffix=".fasta")
#fa
work_dir,temp_dir = split_fasta(input_fa,gene_dict2,root_dir,sub_name="fa")
mv_file(json_data,temp_dir,work_dir,missing_fa,suffix=".fasta")

# bed
parse_gff_to_txt(genome_browser_dir,temp_bed_dir,genome_dict,suffix="gff")
mv_file(json_data,bed_dir,temp_bed_dir,missing_bed,suffix=".bed")