
import csv
import os
import json
import pickle
from tqdm.auto import tqdm
from collections import defaultdict


def read_sam(filepath):
    """
    读取SAM文件并返回行数据。
    """
    rows = []
    with open(filepath, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row[0].startswith("@"):
                rows.append(row)
    return rows


def process_sam_file(filepath, primer_dir):
    sam_data = read_sam(filepath)
    genome_id = os.path.basename(filepath)[:15]
    result = defaultdict(lambda: defaultdict(list))

    for row in sam_data:
        region_name = row[2]  # 区域名
        primer_name = row[0]  # 引物名
        result[region_name][genome_id].append(primer_name)

    json_file = os.path.join(primer_dir, f"{genome_id}.json")
    with open(json_file, "w") as jsonf:
        json.dump(result, jsonf, indent=4)

    return result


def generate_sam_report(sam_dir, output_dir):
    total_report = defaultdict(lambda: defaultdict(list))
    report_dir = os.path.join(output_dir, "primer/nt_sam_info")
    os.makedirs(report_dir, exist_ok=True)
    tsv_file = os.path.join(output_dir, "sam_report.tsv")

    sam_files = sorted(os.listdir(sam_dir))

    with open(tsv_file, "w", newline="") as tsvf:
        writer = csv.writer(tsvf, delimiter="\t")
        writer.writerow(["Genome_id", "Region_name", "Region_num"])

        for sam_file in tqdm(sam_files, total=len(sam_files)):
            sam_path = os.path.join(sam_dir, sam_file)
            genome_id = os.path.basename(sam_file)[:15]
            file_result = process_sam_file(sam_path, report_dir)

            for region_name, genome_entries in file_result.items():
                for genome_id, primers in genome_entries.items():
                    region_num = len(primers)
                    writer.writerow([genome_id, region_name, region_num])
                    total_report[region_name][genome_id].extend(primers)

    # 保存总表为pickle文件
    report_file = os.path.join(output_dir, "bacteria_sam_report.pk")
    with open(report_file, "wb") as pkf:
        pickle.dump(dict(total_report), pkf)

    # 保存总表为JSON文件
    report_file = os.path.join(output_dir, "bacteria_sam_report.json")
    with open(report_file, "w") as jsonf:
        json.dump(total_report, jsonf, indent=4)


if __name__ == "__main__":
    root_dir="/s97/lab/woodman/4w/result"
    sam_dir = f"{root_dir}/merged_sam"

    generate_sam_report(sam_dir, root_dir)
    
    print("Complete!")