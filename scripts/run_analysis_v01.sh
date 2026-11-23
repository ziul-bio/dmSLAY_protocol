#!/usr/bin/bash

set -euo pipefail

bash scripts/merge_reads.sh

bash scripts/trim_reads.sh

python scripts/count_reads.py -i ./data/FASTA_trimmed -o ./results/counts