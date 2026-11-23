#!/usr/bin/bash

#Exit program if any erro accurs or pipe fails
set -euo pipefail


#for i in {1..12}; do
for i in {4,5,6,7,8,9}; do
    echo "Subsampling sample $i"
    seqkit sample -p 0.1 -s 100 data/FASTQ_raw/S${i}_*F.fastq.gz -o data/FASTQ/S${i}_R1.fastq.gz
    seqkit sample -p 0.1 -s 100 data/FASTQ_raw/S${i}_*R.fastq.gz -o data/FASTQ/S${i}_R2.fastq.gz
 
done
