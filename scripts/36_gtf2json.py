import os
import json

def parse_gtf_line(line, file_prefix):
    fields = line.strip().split("\t")
    

    if "gene" not in fields:
        return None
    gene_index = fields.index("gene")

    Region = fields[0].split(" ")[0]
    start = fields[gene_index + 1]
    end = fields[gene_index + 2] 
    gene_id = fields[gene_index + 6].split("gene_id")[1].split('"')[1] 


    return {
        "Gene_id": f"{file_prefix}_{gene_id}",
        "Location": f"{Region}:{start}..{end}"

def convert_gtf_to_jsonl(gtf_file, jsonl_file):
    """Convert a single GTF file to a JSONL file"""
    file_prefix = os.path.basename(gtf_file).replace(".gtf", "")
    with open(gtf_file, "r") as gtf, open(jsonl_file, "w") as jsonl:
        for line in gtf:
            # if line.startswith("#"): 
            #     continue
            record = parse_gtf_line(line, file_prefix)
            if record:
                jsonl.write(json.dumps(record) + "\n")



work_dir = "/s97/lab/woodman/fungi/fungi_mate"
gtf_dir = os.path.join(work_dir,"gtf_dir")
output_file = f"{work_dir}/fungi.gtf.jsonl"

gtf_files = os.listdir(gtf_dir)
with open(output_file, "w") as fungi_jsonl:
    for gtf_file in gtf_files:
        input_path = os.path.join(gtf_dir, gtf_file)
        jsonl_path = os.path.join(work_dir,"gtf2json", gtf_file.replace(".gtf", ".jsonl"))
        

        os.makedirs(os.path.dirname(jsonl_path), exist_ok=True)
        convert_gtf_to_jsonl(input_path, jsonl_path)
        
        with open(jsonl_path, "r") as jsonl:
            fungi_jsonl.writelines(jsonl.readlines())