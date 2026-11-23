#!/usr/bin/bash

#Exit program if any erro accurs or pipe fails
set -euo pipefail

echo "************************************************************************************************"
echo "                          Deep mutational scanning - Part 1 - Count reads"
echo "************************************************************************************************"
echo " "
echo " "

echo "*************************************** Pre-processing reads ***************************************"
echo " "
FASTQ=./data/FASTQ_merged/
OUTDIR1=./data/Fasta
OUTDIR2=./data/FASTA_left_trimmed
OUTDIR3=./data/FASTA_trimmed
OUTDIR4=./results/counts

mkdir -p "$OUTDIR1"
mkdir -p "$OUTDIR2"
mkdir -p "$OUTDIR3"
mkdir -p "$OUTDIR4"



# main loop to process the reads
i=1
while [[ -f "${FASTQ}/S${i}_merged.fastq.gz" ]]; do

    start_time=$(date +%s) # Start time in seconds
    echo "|---------------------- Trimming sample S${i} ----------------------|"

    echo "Filtering reads by amplicon and converting to FASTA"
    seqkit amplicon --threads 50 -m 2 \
        -F CTCCAGCTGCGGGTATCGGAGGAACC \
        -R CATGCCTGCAGGTCGAC \
        "${FASTQ}/S${i}_merged.fastq.gz" | \
        seqkit fq2fa --threads 50 > "${OUTDIR1}/S${i}.fasta"
        echo ""

    # trim adapter from the right and left side of the reads
    # flexbar also has a mismatch option with ratio of 0.1 (10%)
    echo "Trimming reads with flexbar"
    flexbar --threads 50 -ae 0.1 --reads "${OUTDIR1}/S${i}.fasta" --adapters data/adapter_left.fasta --adapter-trim-end LEFT --target "${OUTDIR2}/S${i}"
    flexbar --threads 50 -ae 0.1 --reads "${OUTDIR2}/S${i}.fasta" --adapters data/adapter_right.fasta --adapter-trim-end RIGHT --target "${OUTDIR3}/S${i}"

    end_time=$(date +%s) # End time in seconds
    duration=$(( (end_time - start_time) / 60 ))
    echo "Duration: $duration minutes."
    echo " "

    ((i++))
done

echo "$((i-1)) samples trimmed."

echo "*************************************** Process finished ********************************************"


