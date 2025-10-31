import pysam
import argparse

#print(f"pysam version: {pysam.__version__}")

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file ")
ap.add_argument("-o", "--output", required=True, help="output file")
ap.add_argument("-j", "--outnameinfo", required=True, help="output file of name information")
args = vars(ap.parse_args())
input_fa = args["input"]
output_fq = args["output"]
output_txt = args["outnameinfo"]

# set fragment length
fragment_length = 150

def split_sequence(sequence, fragment_length):
	fragments = [sequence[i:i+fragment_length] for i in range(0, len(sequence), fragment_length)]
	return fragments

with pysam.FastxFile(input_fa) as fh, open(output_fq, "w") as out_fq, open(output_txt, "w") as out_txt:
	for entry in fh:
		sequence_name = entry.name
		sequence = entry.sequence
		fragment_list = split_sequence(sequence, fragment_length)
		name_count = sequence_name + "\t" + str(len(fragment_list)) + "\n"
		out_txt.write(name_count)
		for i in range(0,len(fragment_list)):
			new_name = sequence_name + '/' + str(i)
			new_seq = fragment_list[i]
			new_qual = 'z'*len(new_seq)
			out_fq.write("@%s\n%s\n+\n%s\n" % (new_name, new_seq, new_qual))	