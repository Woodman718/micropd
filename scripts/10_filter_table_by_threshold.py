import numpy as np
import csv, time
from tqdm import tqdm
import os

def try_float(item):
    """Try to convert the input to a float; if it fails, return a very large number."""
    try:
        return float(item)
    except ValueError:
        return float('inf')
    
file_path = '/s1/tangzhidong/project/primerdb/result/PhyloPhlAn/mash_new/Bacteria/bacteria_distmat.tsv'
suf = time.strftime("%m%d_%H%M", time.localtime())
output_path = f'./test/{suf}_bacteria_distmat.tsv'
drop_path = f'./test/{suf}_drop_bacteria.tsv'
threshold_value = 0.001  # Default: 0.001

# Read the first column data from the file
first_column = np.genfromtxt(file_path, delimiter='\t', max_rows=1, dtype='str')[1:]
first_column = np.array([line.split('/')[-1] for line in first_column])

to_remove = set()

with open(file_path, 'r', encoding='utf-8') as file:
    reader = csv.reader(file, delimiter='\t')
    next(reader)  # Skip the header row
    
    for i, row in tqdm(enumerate(reader), desc="Analyzing rows"):
        # Count the number of values that are greater than 0 and less than or equal to the threshold
        count = sum(1 for item in row[1:] if item and try_float(item) <= threshold_value)
        if count >= 2:
            to_remove.add(i)

with open(file_path, 'r', encoding='utf-8') as infile, \
     open(output_path, 'w', newline='', encoding='utf-8') as outfile, \
     open(drop_path, 'w', newline='', encoding='utf-8') as dropfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    drop_writer = csv.writer(dropfile, delimiter='\t')
    
    # Write the header row
    header = next(reader)
    writer.writerow([''] + [first_column[i] for i in range(len(first_column)) if i not in to_remove])
    drop_writer.writerow([''] + [first_column[i] for i in range(len(first_column)) if i in to_remove])
    
    for i, row in tqdm(enumerate(reader), desc="Writing filtered data"):
        if i not in to_remove:
            writer.writerow([first_column[i]] + [row[j + 1] for j in range(len(row) - 1) if j not in to_remove])
        else:
            drop_writer.writerow([first_column[i]] + [row[j + 1] for j in range(len(row) - 1) if j in to_remove])