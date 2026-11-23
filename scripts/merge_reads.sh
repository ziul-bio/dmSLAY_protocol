
#!/usr/bin/env bash

set -euo pipefail

# Usage:
# bash merge_reads_2.sh INPUT_DIR OUTPUT_DIR REPORT_DIR

INPUT_DIR="${1:?Error: provide input directory}"
OUTPUT_DIR="${2:-data/FASTQ_merged}"
REPORT_DIR="${3:-results/fastp_reports}"

# Check dependencies
command -v fastp >/dev/null 2>&1 || {
    echo "Error: fastp is not installed or not in PATH."
    exit 1
}

# Check input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR" "$REPORT_DIR"

# Find and process samples
count=0

shopt -s nullglob

for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do

    # Extract sample name
    filename=$(basename "$R1")
    sample="${filename%_R1.fastq.gz}"

    R2="${INPUT_DIR}/${sample}_R2.fastq.gz"

    # Check matching R2 file
    if [[ ! -f "$R2" ]]; then
        echo "Warning: missing R2 for ${sample}. Skipping."
        continue
    fi

    echo "------------------------------------------------------------"
    echo "Processing sample: ${sample}"
    echo "R1: ${R1}"
    echo "R2: ${R2}"
    echo "------------------------------------------------------------"

    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --merge \
        --out1 "${OUTPUT_DIR}/${sample}_R1_unmerged.fastq.gz" \
        --out2 "${OUTPUT_DIR}/${sample}_R2_unmerged.fastq.gz" \
        --merged_out "${OUTPUT_DIR}/${sample}_merged.fastq.gz" \
        --json "${REPORT_DIR}/${sample}.json" \
        --html "${REPORT_DIR}/${sample}.html"

    ((count+=1))

    echo ""
    echo "============================================================"
    echo ""

done

echo "============================================================"
echo "Processing complete."
echo "Total samples processed: ${count}"
echo "============================================================"
