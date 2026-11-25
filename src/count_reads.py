import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq, BiopythonWarning
from collections import Counter
from pathlib import Path
import hashlib

import warnings
warnings.filterwarnings("ignore", category=BiopythonWarning, message="Partial codon")

import argparse
parser = argparse.ArgumentParser(description='Count the number of reads in a FASTA file')
parser.add_argument('--input_dir', '-i', required=True, help='The directory containing the FASTA files')
parser.add_argument(
    '--output_dir',
    '-o',
    help='The directory to write the results to (default: results/<parent-name>/counts)',
)
parser.add_argument('--minimum-count', '-m', type=int, default=10,
                    help='Minimum total read count to include a peptide (default: 10)')
args = parser.parse_args()

input_dir = Path(args.input_dir)
if not input_dir.is_dir():
    parser.error(f"input directory does not exist: {input_dir}")
output_dir = (
    Path(args.output_dir)
    if args.output_dir
    else Path(__file__).resolve().parent.parent
    / "results"
    / input_dir.parent.name
    / "counts"
)
output_dir.mkdir(parents=True, exist_ok=True)

print('')
print("#####################################################")
print("                    Counting Reads                   ")
print("#####################################################")

# Count nucleotide sequences and aggregate synonymous sequences by peptide.
sample_names = []
peptide_counts = {}

for input_fasta in sorted(input_dir.glob("*.fasta")):
    print(f"Counting reads for {input_fasta.name}.")
    sample = input_fasta.stem
    sample_names.append(sample)

    counts = Counter()
    for record in SeqIO.parse(input_fasta, "fasta"):
        counts[str(record.seq).upper()] += 1

    for sequence, count in counts.items():
        peptide = str(Seq(sequence).translate(to_stop=True))
        peptide_counts.setdefault(peptide, {})[sample] = (
            peptide_counts.setdefault(peptide, {}).get(sample, 0) + count
        )

if not sample_names:
    parser.error(f"no .fasta files found in {input_dir}")

rows = []
for peptide, counts_by_sample in peptide_counts.items():
    row = {"peptide": peptide}
    row.update({sample: counts_by_sample.get(sample, 0) for sample in sample_names})
    rows.append(row)

merged = pl.DataFrame(rows, schema={"peptide": pl.String, **{
    sample: pl.Int64 for sample in sample_names
}})
merged = merged.with_columns(
    pl.sum_horizontal(sample_names).alias("row_sum"),
    pl.col("peptide").str.len_chars().alias("peptide_length"),
)
filtered = merged.filter(pl.col("row_sum") >= args.minimum_count)
filtered = filtered.sort("row_sum", descending=True)


filtered = filtered.with_columns(
    pl.col("peptide")
    .map_elements(
        lambda seq: hashlib.sha256(seq.encode()).hexdigest()[:12],
        return_dtype=pl.String,
    )
    .alias("ID")
)

output_file = output_dir / "peptide_read_counts.txt"
filtered.write_csv(output_file, separator="\t")

print(f"Done! Results saved to {output_file}")
print("\n" * 3, end="")  # adds three blank lines