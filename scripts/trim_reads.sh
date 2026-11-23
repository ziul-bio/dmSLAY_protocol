
#!/usr/bin/env bash

set -euo pipefail

# Usage:
# bash scripts/trim_reads_2.sh \
#     -i data/FASTQ_merged \
#     -o data/FASTA_trimmed \
#     -a data/adapter_left.fasta \
#     -b data/adapter_right.fasta \
#     -t 32 \
#     -f CTCCAGCTGCGGGTATCGGAGGAACC \
#     -r CATGCCTGCAGGTCGAC


# Default parameters
INPUT_DIR=""
TRIMMED_DIR="data/FASTA_trimmed"

ADAPTER_LEFT="data/adapter_left.fasta"
ADAPTER_RIGHT="data/adapter_right.fasta"

THREADS=50

FORWARD_AMPLICON="CTCCAGCTGCGGGTATCGGAGGAACC"
REVERSE_AMPLICON="CATGCCTGCAGGTCGAC"


# Parse arguments
while getopts ":i:o:a:b:t:f:r:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) TRIMMED_DIR="$OPTARG" ;;
        a) ADAPTER_LEFT="$OPTARG" ;;
        b) ADAPTER_RIGHT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f) FORWARD_AMPLICON="$OPTARG" ;;
        r) REVERSE_AMPLICON="$OPTARG" ;;
        h)
            echo "Usage: $0 -i INPUT_DIR [options]"
            echo ""
            echo "Options:"
            echo "  -i  Input FASTQ directory (required)"
            echo "  -o  Output FASTA directory (default: $TRIMMED_DIR)"
            echo "  -a  Left adapter file (default: $ADAPTER_LEFT)"
            echo "  -b  Right adapter file (default: $ADAPTER_RIGHT)"
            echo "  -t  Number of threads (default: $THREADS)"
            echo "  -f  Forward amplicon sequence"
            echo "  -r  Reverse amplicon sequence"
            echo "  -h  Display this help message"
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument."
            exit 1
            ;;
        \?)
            echo "Error: invalid option -$OPTARG"
            exit 1
            ;;
    esac
done


# Validate required arguments
if [[ -z "$INPUT_DIR" ]]; then
    echo "Error: input directory (-i) is required."
    exit 1
fi


# Validate input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory does not exist: $INPUT_DIR"
    exit 1
fi


# Validate thread count
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer."
    exit 1
fi


# Check dependencies
for cmd in seqkit flexbar; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    }
done


# Check adapter files
for adapter in "$ADAPTER_LEFT" "$ADAPTER_RIGHT"; do
    if [[ ! -f "$adapter" ]]; then
        echo "Error: adapter file not found: $adapter"
        exit 1
    fi
done


# Create output directory
mkdir -p "$TRIMMED_DIR"


# Create temporary directory and configure cleanup
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT


# Identify input files
shopt -s nullglob
FASTQ_FILES=("$INPUT_DIR"/*_merged.fastq.gz)

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "Error: no *_merged.fastq.gz files found in $INPUT_DIR"
    exit 1
fi


# Pipeline information
echo "************************************************************"
echo "              Deep Mutational Scanning"
echo "                    Read Trimming"
echo "************************************************************"
echo "Input directory: $INPUT_DIR"
echo "Temporary directory: $TEMP_DIR"
echo "Output directory: $TRIMMED_DIR"
echo "Threads: $THREADS"
echo "Samples found: ${#FASTQ_FILES[@]}"
echo ""


# Process samples
count=0

for FASTQ in "${FASTQ_FILES[@]}"; do

    start_time=$(date +%s)

    # Extract sample name
    filename=$(basename "$FASTQ")
    sample="${filename%_merged.fastq.gz}"

    echo "------------------------------------------------------------"
    echo "Processing sample: $sample"
    echo "------------------------------------------------------------"


    # Step 1: Filter reads by amplicon and convert to FASTA
    echo "Filtering reads by amplicon and converting to FASTA..."

    seqkit amplicon --threads "$THREADS" -m 2 \
        -F "$FORWARD_AMPLICON" \
        -R "$REVERSE_AMPLICON" \
        "$FASTQ" | \
        seqkit fq2fa --threads "$THREADS" \
        > "${TEMP_DIR}/${sample}.fasta"


    # Step 2: Trim left adapter
    echo "Trimming left adapter..."
    flexbar --threads "$THREADS" -ae 0.1 \
        --reads "${TEMP_DIR}/${sample}.fasta" \
        --adapters "$ADAPTER_LEFT" \
        --adapter-trim-end LEFT \
        --target "${TEMP_DIR}/${sample}_left_tmp"


    # Validate intermediate output
    if [[ ! -f "${TEMP_DIR}/${sample}_left_tmp.fasta" ]]; then
        echo "Error: Flexbar did not generate the left-trimmed FASTA for $sample."
        exit 1
    fi


    # Step 3: Trim right adapter
    echo "Trimming right adapter..."
    flexbar --threads "$THREADS" -ae 0.1 \
        --reads "${TEMP_DIR}/${sample}_left_tmp.fasta" \
        --adapters "$ADAPTER_RIGHT" \
        --adapter-trim-end RIGHT \
        --target "${TRIMMED_DIR}/${sample}"


    # Processing time
    end_time=$(date +%s)
    duration=$(( (end_time - start_time) / 60 ))

    ((count+=1))

    echo "Finished processing $sample."
    echo "Duration: $duration minutes."
    echo ""

done


# Final report
echo "************************************************************"
echo "Trimming complete."
echo "Total samples processed: $count"
echo "************************************************************"
