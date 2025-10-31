import os
import argparse
import pysam
from tqdm import tqdm
# from collections import defaultdict

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Split merged fasta files into genome-specific fasta files.")
parser.add_argument("-m", "--merged_fasta_dir", required=True, help="Directory of merged fasta files")
parser.add_argument("-g", "--gene_name_dir", required=True, help="Directory of gene name txt files")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory for genome-specific fasta files")
args = parser.parse_args()

merged_fasta_dir = args.merged_fasta_dir
gene_name_dir = args.gene_name_dir
output_dir = args.output_dir

os.makedirs(output_dir, exist_ok=True)

# Load merged fasta files into a dictionary, using pysam for efficient reading
def load_fasta_to_dict(merged_fasta_dir):
    fasta_dict = {}
    print("Loading merged fasta files into memory...")
    
    # Loop through all merged fasta files
    for fasta_file in os.listdir(merged_fasta_dir):
        fasta_path = os.path.join(merged_fasta_dir, fasta_file)
        
        with pysam.FastxFile(fasta_path, "r") as fi_fa:
            for entry in fi_fa:
                fasta_dict[entry.name] = entry.sequence
    
    print(f"Loaded {len(fasta_dict)} genes into memory.")
    return fasta_dict

# Process each gene ID file and generate genome-specific fasta files
def process_gene_name_files(gene_file, fasta_dict, fasta_ids, output_dir):
    genome_name = gene_file.replace('.txt', '')
    output_fasta = os.path.join(output_dir, f"{genome_name}.fasta")

    with open(os.path.join(gene_name_dir, gene_file), 'r') as gf, open(output_fasta, 'w') as of:
        gene_ids = {line.strip() for line in gf} # Read gene names into a set

        # Find the intersection of gene IDs with the fasta dictionary keys and sort them
        matched_ids = sorted(gene_ids.intersection(fasta_ids))
        
        for gene_id in matched_ids:
            of.write(f">{gene_id}\n")
            of.write(f"{fasta_dict[gene_id]}\n")
    
    # print(f"Written {len(matched_ids)} genes to {output_fasta}.")

# Main function
def main():
    # Load merged fasta files into a dictionary
    fasta_dict = load_fasta_to_dict(merged_fasta_dir)
    fasta_ids = set(fasta_dict.keys())
    
    gene_files = os.listdir(gene_name_dir)
    
    # Process gene ID files and generate new fasta files
    for gene_file in tqdm(gene_files, desc="Processing"):
        process_gene_name_files(gene_file, fasta_dict, fasta_ids, output_dir)

if __name__ == "__main__":
    main()