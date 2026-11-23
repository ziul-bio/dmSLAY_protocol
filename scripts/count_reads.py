import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

import argparse
parser = argparse.ArgumentParser(description='Count the number of reads in a FASTA file')
parser.add_argument('--input_dir', '-i', required=True, help='The directory containing the FASTA files')
parser.add_argument('--output_dir','-o', required=True, help='The directory to write the results to')
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

# Step 1: Read the FASTA file
sample_counts = []

for file in os.listdir(args.input_dir):
    if file.endswith(".fasta"):
        print(f"Counting reads for {file}.")

        input_fasta = os.path.join(args.input_dir, file)
        sample = file.removesuffix(".fasta")

        # Count sequences
        counts = Counter()
        for record in SeqIO.parse(input_fasta, "fasta"):
            counts[str(record.seq)] += 1

        # Store counts as a Series
        sample_counts.append(pd.Series(counts, name=sample))

# Combine all samples
merged = pd.concat(sample_counts, axis=1).fillna(0).astype(int)
merged = merged.reset_index()
merged = merged.rename(columns={'index': 'Sequence'})
print(merged.head())

# Translate sequences to peptides and group by peptide
merged['peptide'] = merged['Sequence'].apply(lambda x: "".join(Seq(x).translate(to_stop=True)))
merged = merged.groupby(["peptide"], as_index=False).sum().copy()

# Compute row sums and peptide lengths
merged['row_sum'] = merged.iloc[:, 2:].sum(axis=1)
merged.sort_values(by='row_sum', ascending=False, inplace=True)
merged['peptide_length'] = merged['peptide'].apply(len)


# Save results, filtering out peptides with low read counts
merged.query('row_sum > 10').to_csv(f"{args.output_dir}/peptide_read_counts.txt", sep="\t")