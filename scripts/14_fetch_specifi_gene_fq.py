import pysam
import argparse

# Define the argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file from fasta")
ap.add_argument("-r", "--background", required=True, help="input txt file include gene names")
ap.add_argument("-o", "--output", required=True, help="output file for sequences that match the criteria")
ap.add_argument("-f", "--filter", required=True, help="output file for sequences that do not match the criteria")

# Parse the arguments
args = vars(ap.parse_args())

# Assign arguments to variables
input_fa = args["input"]
input_txt = args["background"]
output_fa = args["output"]
filter_fa = args["filter"]

# Read the list of gene IDs
with open(input_txt, 'r') as fi_txt:
    gene_ids = [line.strip() for line in fi_txt]

# Process the FASTA file
with pysam.FastxFile(input_fa, "r") as fi_fa, open(output_fa, "w") as fo, open(filter_fa, "w") as ff:
    for entry in fi_fa:
        sequence_name = entry.name
        sequence = entry.sequence
        if sequence_name not in gene_ids:
            # Write sequences that do not match the criteria to the filter file
            ff.write(">" + sequence_name + "\n")
            ff.write(sequence + "\n")
        else:
            # Write sequences that match the criteria to the output file
            fo.write(">" + sequence_name + "\n")
            fo.write(sequence + "\n")