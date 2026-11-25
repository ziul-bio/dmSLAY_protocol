#!/usr/bin/bash

#Exit program if any erro accurs or pipe fails
set -euo pipefail


#for i in {1..12}; do
for i in {1,2,3,10,11,12}; do
    echo "Subsampling sample $i"
    seqkit sample -p 0.05 -s 100 data/FASTQ_raw/S${i}_*F.fastq.gz -o data/FASTQ_example/S${i}_R1.fastq.gz
    seqkit sample -p 0.05 -s 100 data/FASTQ_raw/S${i}_*R.fastq.gz -o data/FASTQ_example/S${i}_R2.fastq.gz
 
done
