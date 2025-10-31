import os
import pandas as pd
from Bio.Seq import Seq


if len(sys.argv) != 5:
    print("Usage: python filename.py <input_file> <output_dir> <intermediate_dir> <csv_name>")
    sys.exit(1)
    
data_path = sys.argv[1]
outpu_dir = sys.argv[2]
intermediate_dir = sys.argv[3]
csv_name = sys.argv[4]

out_path = f"{outpu_dir}/{csv_name}_score.csv"
out_intermediate = f"{out_intermediate}/{csv_name}_intermediate.csv"

# Function to calculate Gibbs free energy for a sequence
def calculate_gibbs_energy(seq):
    """
    Calculate Gibbs free energy for a given sequence.
    
    seq: The DNA sequence (string).
    Returns the calculated Gibbs energy (float).
    """
    seq = Seq(seq)
    # Calculate 3' end Gibbs energy
    subseq_3prime = str(seq[-6:])  # Last 6 bases
    score_3prime = 0
    for i in range(len(subseq_3prime) - 1):
        dimer = "/".join([subseq_3prime[i:i+2], str(Seq(subseq_3prime[i:i+2]).complement())])
        score_3prime += deltaG_dict.get(dimer, 0)
    # Adjust score based on the last base
    sixth_base = subseq_3prime[5]
    score_3prime += 0.98 if sixth_base in ("A", "T") else 1.00

    # Calculate 5' end Gibbs energy
    subseq_5prime = str(seq[:6])  # First 6 bases
    score_5prime = 0
    for i in range(len(subseq_5prime) - 1):
        dimer = "/".join([subseq_5prime[i:i+2], str(Seq(subseq_5prime[i:i+2]).complement())])
        score_5prime += deltaG_dict.get(dimer, 0)
    # Adjust score based on the first base
    first_base = subseq_5prime[0]
    score_5prime += 0.98 if first_base in ("A", "T") else 1.00

    # Total Gibbs free energy
    gibbs_energy = score_3prime - score_5prime
    return gibbs_energy

# Calculate GC content
def calculate_gc_content(sequence):
    """
    Calculate GC content percentage.
    
    sequence: DNA sequence (string).
    Returns the GC content as a percentage (float).
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) * 100


# Load data
data = pd.read_csv(data_path)
data['KINGDOM'] = data['KINGDOM'].fillna("NA")

deltaG_dict = {
    "AA/TT": -0.73,
    "AT/TA": -0.61,
    "AC/TG": -1.16,
    "AG/TC": -0.92,
    "TA/AT": -0.32,
    "TT/AA": -0.73,
    "TC/AG": -1.03,
    "TG/AC": -1.16,
    "CA/GT": -1.16,
    "CT/GA": -0.92,
    "CC/GG": -1.57,
    "CG/GC": -1.81,
    "GA/CT": -1.03,
    "GT/CA": -1.16,
    "GC/CG": -1.92,
    "GG/CC": -1.57
}

# Initialize intermediate variables
results = []  # Store rows with calculated scores and intermediate values

# Process each row
for _, row in data.iterrows():
    # Calculate Gibbs energy for left and right primers
    gibbs_left = calculate_gibbs_energy(row['PRIMER_LEFT_SEQUENCE'])
    gibbs_right = calculate_gibbs_energy(row['PRIMER_RIGHT_SEQUENCE'])
    total_gibbs = gibbs_left + gibbs_right

    # Calculate GC content and length
    gc_left = calculate_gc_content(row['PRIMER_LEFT_SEQUENCE'])
    gc_right = calculate_gc_content(row['PRIMER_RIGHT_SEQUENCE'])
    len_left = len(row['PRIMER_LEFT_SEQUENCE'])
    len_right = len(row['PRIMER_RIGHT_SEQUENCE'])

    # Calculate melting temperatures (Tm)
    tm_left = 59.9 + (0.41 * gc_left) - (675 / len_left)
    tm_right = 59.9 + (0.41 * gc_right) - (675 / len_right)
    delta_tm = abs(tm_left - tm_right)

    # Calculate hairpin and self-complementarity
    hairpin = row['PRIMER_LEFT_HAIRPIN_TH'] + row['PRIMER_RIGHT_HAIRPIN_TH']
    self_comp = (
        row['PRIMER_LEFT_SELF_ANY_TH'] +
        row['PRIMER_RIGHT_SELF_ANY_TH'] +
        row['PRIMER_PAIR_COMPL_END_TH']
    )

    # Calculate final score
    score = delta_tm + total_gibbs + hairpin + self_comp

    # Append results to intermediate table
    results.append({
        'GENOME_ID': row['GENOME_ID'],
        'GENE_ID': row['GENE_ID'],
        'PRIMER_LEFT_SEQUENCE': row['PRIMER_LEFT_SEQUENCE'],
        'PRIMER_RIGHT_SEQUENCE': row['PRIMER_RIGHT_SEQUENCE'],
        'score': score,
        'gc_contents_L': gc_left,
        'gc_contents_R': gc_right,
        'Tm_L': tm_left,
        'Tm_R': tm_right,
        'delta_Tm': delta_tm,
        'hairpin': hairpin,
        'self': self_comp,
        'G': total_gibbs
    })

intermediate_df = pd.DataFrame(results)
data['score'] = intermediate_df['score']


intermediate_df.to_csv(out_intermediate, index=False)
data.to_csv(out_path, index=False)