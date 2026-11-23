#!/usr/bin/bash

set -euo pipefail

# bash scripts/merge_reads_2.sh ./data/FASTQ

# bash scripts/trim_reads_2.sh \
#     -i data/FASTQ_merged \
#     -o data/FASTA_trimmed \
#     -a data/adapter_left.fasta \
#     -b data/adapter_right.fasta \
#     -t 32 \
#     -f CTCCAGCTGCGGGTATCGGAGGAACC \
#     -r CATGCCTGCAGGTCGAC

python scripts/count_reads_2.py -i ./data/FASTA_trimmed -o ./results/counts2