import os
import pandas as pd
from Bio import SeqIO
from collections import Counter

import argparse
parser = argparse.ArgumentParser(description='Count the number of reads in a FASTA file')
parser.add_argument('--input_dir', '-i', required=True, help='The directory containing the FASTA files')
parser.add_argument('--output_dir','-o', required=True, help='The directory to write the results to')
args = parser.parse_args()


# Step 1: Read the FASTA file
input_dir = args.input_dir
output_dir = args.output_dir

for file in os.listdir(input_dir):
    if file.endswith(".fasta"):
        print(f"Counting reads for {file}.")
        input_fasta = os.path.join(input_dir, file)

        ## Step 2: Count sequences
        counts = Counter()
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            counts[seq] += 1

        # Step 3: Output the results
        sample = file.strip('.fasta')
        print(f'Writing results to sample: {sample}.')
        out_path = f"{output_dir}/{sample}.txt"
    
        with open(f"{out_path}", "w") as output_file:
            output_file.write(f"Sequence\t{sample}\n")
            for sequence, count in counts.items():
                output_file.write(f"{sequence}\t{count}\n")
