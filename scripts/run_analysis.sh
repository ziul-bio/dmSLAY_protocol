#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT_DIR=""
MERGED_DIR=""
TRIMMED_DIR=""
RESULTS_DIR=""
REPORT_DIR=""
ADAPTER_LEFT=""
ADAPTER_RIGHT=""
FIVE_PRIME_FLANK=""
THREE_PRIME_FLANK=""
THREADS=32
MINIMUM_COUNT=10
COLDATA_FILE=""
DESEQ2_REFERENCE=""
DESEQ2_TREATMENT=""
DESEQ2_OUTPUT_DIR=""
DESEQ2_MINIMUM_COUNT=30

usage() {
    cat <<EOF
Usage: bash scripts/run_analysis.sh --input-dir DIR --five-prime-sequence SEQ \\
  --three-prime-sequence SEQ --coldata-file FILE [options]

Required arguments:
  -i, --input-dir DIR              Directory containing paired reads
  -5, --five-prime-sequence SEQ    5-prime/forward flanking sequence
  -3, --three-prime-sequence SEQ   3-prime/reverse flanking sequence
  --coldata-file FILE              TSV/CSV with sample and group/condition columns

Options:
  --merged-dir DIR                 Merged FASTQ output directory
  --trimmed-dir DIR                Trimmed FASTA output directory
  --results-dir DIR                Count output directory
  --report-dir DIR                 fastp report directory
  -a, --adapter-left FILE          Optional left adapter FASTA override
  -b, --adapter-right FILE         Optional right adapter FASTA override
  -f, --forward-sequence SEQ   Alias for 5-prime flanking sequence
  -r, --reverse-sequence SEQ   Alias for 3-prime flanking sequence
  -t, --threads N                  Worker threads
  -m, --minimum-count N            Minimum peptide read count
  --deseq2-reference GROUP         DESeq2 reference group
  --deseq2-treatment GROUP         DESeq2 treatment group(s), comma-separated
  --deseq2-output-dir DIR          DESeq2 output directory
  --deseq2-minimum-count N         Minimum peptide count for DESeq2 (default: 30)
  -h, --help                       Display this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir) INPUT_DIR="$2"; shift 2 ;;
        --merged-dir) MERGED_DIR="$2"; shift 2 ;;
        --trimmed-dir) TRIMMED_DIR="$2"; shift 2 ;;
        --results-dir) RESULTS_DIR="$2"; shift 2 ;;
        --report-dir) REPORT_DIR="$2"; shift 2 ;;
        -a|--adapter-left) ADAPTER_LEFT="$2"; shift 2 ;;
        -b|--adapter-right) ADAPTER_RIGHT="$2"; shift 2 ;;
        -5|--five-prime-sequence) FIVE_PRIME_FLANK="$2"; shift 2 ;;
        -3|--three-prime-sequence) THREE_PRIME_FLANK="$2"; shift 2 ;;
        -f|--forward-sequence) FIVE_PRIME_FLANK="$2"; shift 2 ;;
        -r|--reverse-sequence) THREE_PRIME_FLANK="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -m|--minimum-count) MINIMUM_COUNT="$2"; shift 2 ;;
        --coldata-file|--coldata) COLDATA_FILE="$2"; shift 2 ;;
        --deseq2-reference|--reference) DESEQ2_REFERENCE="$2"; shift 2 ;;
        --deseq2-treatment|--treatment) DESEQ2_TREATMENT="$2"; shift 2 ;;
        --deseq2-output-dir) DESEQ2_OUTPUT_DIR="$2"; shift 2 ;;
        --deseq2-minimum-count) DESEQ2_MINIMUM_COUNT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

missing_args=()
[[ -n "$INPUT_DIR" ]] || missing_args+=(--input-dir)
[[ -n "$FIVE_PRIME_FLANK" ]] || missing_args+=(--five-prime-sequence)
[[ -n "$THREE_PRIME_FLANK" ]] || missing_args+=(--three-prime-sequence)
[[ -n "$COLDATA_FILE" ]] || missing_args+=(--coldata-file)
if ((${#missing_args[@]} > 0)); then
    printf 'Error: missing required arguments: %s\n' "${missing_args[*]}" >&2
    printf 'Run with --help to see usage.\n' >&2
    exit 2
fi

INPUT_PARENT_DIR="$(dirname "$INPUT_DIR")"
MERGED_DIR="${MERGED_DIR:-$INPUT_PARENT_DIR/FASTQ_merged}"
TRIMMED_DIR="${TRIMMED_DIR:-$INPUT_PARENT_DIR/FASTA_trimmed}"
RESULTS_DIR="${RESULTS_DIR:-$INPUT_PARENT_DIR/results/counts}"
REPORT_DIR="${REPORT_DIR:-$INPUT_PARENT_DIR/results/fastp_reports}"
DESEQ2_OUTPUT_DIR="${DESEQ2_OUTPUT_DIR:-$INPUT_PARENT_DIR/results/deseq2}"

python3 "$SCRIPT_DIR/merge_reads.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$MERGED_DIR" \
    --report-dir "$REPORT_DIR"

bash "$SCRIPT_DIR/trim_reads.sh" \
    -i "$MERGED_DIR" \
    -o "$TRIMMED_DIR" \
    -a "$ADAPTER_LEFT" \
    -b "$ADAPTER_RIGHT" \
    -t "$THREADS" \
    -5 "$FIVE_PRIME_FLANK" \
    -3 "$THREE_PRIME_FLANK"

python3 "$SCRIPT_DIR/count_reads.py" \
    --input_dir "$TRIMMED_DIR" \
    --output_dir "$RESULTS_DIR" \
    --minimum-count "$MINIMUM_COUNT"

DESEQ2_ARGS=(
    "$SCRIPT_DIR/deseq2.r"
    --count-file "$RESULTS_DIR/peptide_read_counts.txt"
    --output-dir "$DESEQ2_OUTPUT_DIR"
    --minimum-count "$DESEQ2_MINIMUM_COUNT"
    --coldata-file "$COLDATA_FILE"
)
if [[ -n "$DESEQ2_REFERENCE" ]]; then
    DESEQ2_ARGS+=(--reference "$DESEQ2_REFERENCE")
fi
if [[ -n "$DESEQ2_TREATMENT" ]]; then
    DESEQ2_ARGS+=(--treatment "$DESEQ2_TREATMENT")
fi
Rscript "${DESEQ2_ARGS[@]}"