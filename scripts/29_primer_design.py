import os
import argparse
import pysam
import primer3
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input fasta file of specific genes")
ap.add_argument("-o", "--output", required=True, help="path of primer output file")
ap.add_argument("--log", help="path of primer log file")
args = vars(ap.parse_args())
input_fa = args["input"]
output_dir = args["output"]
log_dir = args["log"]


def primer_design(sequence_id, sequence_template, filename):
	# sequence and argument
	seq_args = {
		'SEQUENCE_ID': sequence_id,
		'SEQUENCE_TEMPLATE': sequence_template,
		'SEQUENCE_INCLUDED_REGION': [0,len(sequence_template)]
	}
	# primer condition
	global_args = {
		'PRIMER_NUM_RETURN': 5,
		'PRIMER_OPT_SIZE': 21,
		'PRIMER_MIN_SIZE': 18,
		'PRIMER_MAX_SIZE': 25,
		'PRIMER_OPT_TM': 59.0,
		'PRIMER_MIN_TM': 57.0,
		'PRIMER_MAX_TM': 63.0,
		'PRIMER_MIN_GC': 30.0,
		'PRIMER_MAX_GC': 60.0,
		'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
		'PRIMER_MAX_POLY_X': 100,
		'PRIMER_INTERNAL_MAX_POLY_X': 100,
		'PRIMER_SALT_MONOVALENT': 50.0,
		'PRIMER_DNA_CONC': 50.0,
		'PRIMER_MAX_NS_ACCEPTED': 0,
		'PRIMER_MAX_SELF_ANY': 8,
		'PRIMER_MAX_SELF_END': 3,
		'PRIMER_PAIR_MAX_COMPL_ANY': 8,
		'PRIMER_PAIR_MAX_COMPL_END': 3,
		'PRIMER_PRODUCT_SIZE_RANGE': [100,500],
		'PRIMER_GC_CLAMP': 1
	}

	primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
	primer3_result_table_dict = {}
	for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
		primer_id = str(i)
		for key in primer3_result:
			if primer_id in key:
				info_tag = key.replace("_" + primer_id, "")
				try:
					primer3_result_table_dict[info_tag]
				except:
					primer3_result_table_dict[info_tag] = []
				finally:
					primer3_result_table_dict[info_tag].append(primer3_result[key])

	index = []

	for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
		index.append("PRIMER_PAIR_" + str(i))

	primer3_result_df = pd.DataFrame(primer3_result_table_dict, index=index)
	primer3_result_df = primer3_result_df.T
	primer3_result_df.to_csv(filename + "_primer3_result.csv", sep=",")


genome_id = os.path.basename(input_fa)[:15]
prefix = os.path.join(output_dir,genome_id)
prelog = os.path.join(log_dir,genome_id+".log")

with pysam.FastxFile(input_fa) as fh:
	for entry in fh:
		sequence_id = entry.name
		sequence = entry.sequence
		csv_name = prefix + "_" + sequence_id
  
		if len(sequence) < 100:
			print(f"Sequence length for {sequence_id} is too short: {len(sequence)}")
			with open(prelog, 'a') as f:
				f.write(f"Sequence length for {sequence_id} is too short: {len(sequence)}\n")
			continue

		primer_design(sequence_id, sequence, csv_name)