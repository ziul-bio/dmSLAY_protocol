#!/usr/bin/bash

# merge reads
INPUT_DIR="data/FASTQ"
OUTPUT_DIR="data/FASTQ_merged"
REPORT_DIR="results/fastp_reports"
mkdir -p "$OUTPUT_DIR" "$REPORT_DIR"

i=1
while [[ -f "${INPUT_DIR}/S${i}_R1.fastq.gz" && -f "${INPUT_DIR}/S${i}_R2.fastq.gz" ]]; do
    
    R1="${INPUT_DIR}/S${i}_R1.fastq.gz"
    R2="${INPUT_DIR}/S${i}_R2.fastq.gz"

    echo "|---------------------- Merging sample S${i} ----------------------|"
    echo "R1: ${R1} --> R2: ${R2}"
    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --merge \
        --out1 "${OUTPUT_DIR}/S${i}_R1_unmerged.fastq.gz" \
        --out2 "${OUTPUT_DIR}/S${i}_R2_unmerged.fastq.gz" \
        --merged_out "${OUTPUT_DIR}/S${i}_merged.fastq.gz" \
        --json "${REPORT_DIR}/S${i}.json" \
        --html "${REPORT_DIR}/S${i}.html"
    echo "|-------------------------------END-------------------------------|"
    echo ""
    echo ""
    
    ((i++))
done

echo "$((i-1)) samples processed."