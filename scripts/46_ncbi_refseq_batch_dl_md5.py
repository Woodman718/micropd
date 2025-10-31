import subprocess
import os
import time
from concurrent.futures import ThreadPoolExecutor
import tqdm.auto as tqdm

root_dir = "/home/caisongbo/s3/json"
# id_file = os.path.join(root_dir, "speciGenome_38084.txt")
id_file = os.path.join(root_dir, "failed_ids2.txt")
zip_dir = os.path.join(root_dir, "ncbi_bacteria_refseq")
jsonl_dir = os.path.join(root_dir, "sequence_reports")
md5_dir = os.path.join(root_dir, "MD5-3")
failed_md5_file = os.path.join(root_dir, "failed_md5.txt")
failed_file = os.path.join(root_dir, "failed_ids3.txt")

os.makedirs(zip_dir, exist_ok=True)
os.makedirs(jsonl_dir, exist_ok=True)
os.makedirs(md5_dir, exist_ok=True)

with open(id_file, 'r') as f:
    ids = [line.strip() for line in f]

failed_ids = []
md5_count = 1
nw = 10

def download_and_validate(acc_id, current_batch):
    zip_file = os.path.join(zip_dir, f"{acc_id}.zip")
    jsonl_file = os.path.join(jsonl_dir, f"{acc_id}_seq_report.jsonl")
    temp_md5_file = os.path.join(zip_dir, f"{acc_id}_md5sum.txt")
    
    if os.path.exists(jsonl_file):
        return f"{acc_id} already processed. Skipping..."
    
    command = f"datasets download genome accession {acc_id} --assembly-source RefSeq --include seq-report --no-progressbar --filename {zip_file}"
    try:
        # 下载文件
        subprocess.run(command, shell=True, check=True)
        
        # 解压缩 md5sum.txt 和 sequence_report.jsonl
        extract_md5_command = f"unzip -p {zip_file} md5sum.txt > {temp_md5_file}"
        subprocess.run(extract_md5_command, shell=True, check=True)
        
        extract_jsonl_command = f"unzip -p {zip_file} ncbi_dataset/data/{acc_id}/sequence_report.jsonl > {jsonl_file}"
        subprocess.run(extract_jsonl_command, shell=True, check=True)
        
        # 找到sequence_report.jsonl对应的MD5
        expected_md5 = None
        with open(temp_md5_file, 'r') as f:
            for line in f:
                filename = line.split()[1]
                if "sequence_report.jsonl" in filename:
                    expected_md5 = line.split()[0]
                    break
        
        if expected_md5 is None:
            raise ValueError(f"MD5 for sequence_report.jsonl not found in {temp_md5_file}")
        
        # 将MD5值记录到批次文件
        md5_file_path = os.path.join(md5_dir, f"md5.{current_batch}.txt")
        with open(md5_file_path, 'a') as md5_file:
            md5_file.write(f"{expected_md5}  {jsonl_file}\n")
        
        # 清理临时文件
        os.remove(temp_md5_file)
        return f"{acc_id} processed and verified successfully."
    except (subprocess.CalledProcessError, ValueError) as e:
        failed_ids.append(acc_id)
        return f"{acc_id} failed: {e}"

def validate_md5_for_batch(batch_ids, current_batch):
    md5_file_path = os.path.join(md5_dir, f"md5.{current_batch}.txt")
    failed_batch_ids = []
    try:
        check_md5_command = f"md5sum -c {md5_file_path} --quiet"
        subprocess.run(check_md5_command, shell=True, check=True)
        return f"Batch {current_batch} MD5 verification passed."
    except subprocess.CalledProcessError as e:
        # 如果有验证失败，找出哪些文件失败并记录
        with open(md5_file_path, 'r') as md5_file:
            for line in md5_file:
                file_name, status = line.split()
                if status != 'OK':
                    failed_batch_ids.append(file_name)

        # 记录验证失败的文件
        if failed_batch_ids:
            with open(failed_md5_file, "a") as f:
                for failed_id in failed_batch_ids:
                    f.write(f"{failed_id} MD5 verification failed\n")
        return f"Batch {current_batch} MD5 verification failed for {len(failed_batch_ids)} files."


# 并发下载和处理
print("Starting downloads and verification...")
with ThreadPoolExecutor(max_workers=nw) as executor:
    for i in tqdm.tqdm(range(0, len(ids), 10), desc="Processing Batches"):
        batch_ids = ids[i:i+10]
        
        # 批次处理
        results = list(executor.map(download_and_validate, batch_ids, [md5_count] * len(batch_ids)))
        md5_count += 1 
        
        # 验证当前批次的 MD5 校验
        md5_verification_result = validate_md5_for_batch(batch_ids, md5_count - 1)
        print(md5_verification_result)

# 打印失败的 ID 列表
if failed_ids:
    with open(failed_file, "w") as f:
        f.writelines("\n".join(failed_ids))
    print(f"Failed IDs saved to {failed_file}")

print("Download and processing completed.")