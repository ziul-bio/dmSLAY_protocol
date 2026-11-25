
#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Usage:
# bash src/trim_reads.sh \
#     -i data/FASTQ_merged \
#     -o data/FASTA_trimmed \
#     -t 8 \
#     -5 CTCCAGCTGCGGGTATCGGAGGAACC \
#     -3 GTCGACCTGCAGGCATG


# Default parameters
THREADS=4
INPUT_DIR=""
TRIMMED_DIR=""

FIVE_PRIME_FLANK=""
THREE_PRIME_FLANK=""

# Parse arguments
while getopts ":i:o:t:f:r:5:3:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) TRIMMED_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f|5) FIVE_PRIME_FLANK="$OPTARG" ;;
        r|3) THREE_PRIME_FLANK="$OPTARG" ;;
        h)
            echo "Usage: $0 -i INPUT_DIR [options]"
            echo ""
            echo "Options:"
            echo "  -i  Input FASTQ directory (required)"
            echo "  -o  Output FASTA directory (default: results/<parent-name>/FASTA_trimmed)"
            echo "  -t  Number of threads (default: $THREADS)"
            echo "  -5  5-prime (forward) flanking sequence (required)"
            echo "  -3  3-prime (reverse) flanking sequence (required)"
            echo "  -f  Alias for -5"
            echo "  -r  Alias for -3"
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

normalize_sequence() {
    printf '%s' "$1" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]'
}

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

FIVE_PRIME_FLANK="$(normalize_sequence "$FIVE_PRIME_FLANK")"
THREE_PRIME_FLANK="$(normalize_sequence "$THREE_PRIME_FLANK")"
if [[ -z "$FIVE_PRIME_FLANK" || -z "$THREE_PRIME_FLANK" ]]; then
    echo "Error: both 5-prime (-5/-f) and 3-prime (-3/-r) sequences are required." >&2
    exit 1
fi

if [[ -z "$TRIMMED_DIR" ]]; then
    INPUT_NAME="$(basename "$(dirname "$INPUT_DIR")")"
    TRIMMED_DIR="$PROJECT_DIR/results/$INPUT_NAME/FASTA_trimmed"
fi


# Validate thread count
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer."
    exit 1
fi


# Check dependencies
for cmd in seqkit; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    }
done


# Create output directory
mkdir -p "$TRIMMED_DIR"


reverse_complement() {
    python3 - "$1" <<'PY'
import sys

def reverse_complement(sequence):
    """
    # Calculate the reverse complement of a DNA/RNA sequence.

    # Supports standard nucleotides (A, C, G, T, U) and IUPAC
    # nucleotide ambiguity codes:
    #                     R = A or G (purine)
    #                     Y = C or T (pyrimidine)
    #                     K = G or T
    #                     M = A or C
    #                     S = G or C
    #                     W = A or T
    #                     B = C, G, or T
    #                     V = A, C, or G
    #                     D = A, G, or T
    #                     H = A, C, or T
    #                     N = any nucleotide
    # Parameters
    # ----------
    # sequence : str
    #     Input nucleotide sequence.

    # Returns
    # -------
    # str
    #     Reverse complement of the input sequence.

    # Raises
    # ------
    # ValueError
    #     If the sequence contains invalid nucleotide characters.
    """

    sequence = sequence.upper()

    valid_bases = set("ACGTURYKMSWBVDHN")

    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(
            f"Invalid nucleotide characters: {', '.join(sorted(invalid_bases))}"
        )

    complement = str.maketrans(
        "ACGTURYKMSWBVDHN",
        "TGCAAYRMKSWVB DHN".replace(" ", ""),
    )

    return sequence.translate(complement)[::-1]


try:
    print(reverse_complement(sys.argv[1]))
except ValueError as e:
    sys.exit(f"Error: {e}")

PY
}

# Compute the reverse complement of the 3' flank.
REVERSE_AMPLICON="$(reverse_complement "$THREE_PRIME_FLANK")"

# Identify input files
shopt -s nullglob
FASTQ_FILES=()
for pattern in "*_merged.fastq.gz" "*_merged.fq.gz" "*_merged.fastq" "*_merged.fq"; do
    FASTQ_FILES+=("$INPUT_DIR"/$pattern)
done

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "Error: no merged FASTQ files found in $INPUT_DIR"
    exit 1
fi


# Pipeline information
echo "************************************************************"
echo "                 Amplicon Extraction"
echo "************************************************************"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $TRIMMED_DIR"
echo "Threads: $THREADS"
echo "Samples found: ${#FASTQ_FILES[@]}"
echo "5-prime flank: $FIVE_PRIME_FLANK"
echo "3-prime flank: $THREE_PRIME_FLANK"
echo "3-prime reverse complement for seqkit -R: $REVERSE_AMPLICON"
echo ""


# Process samples
count=0

for FASTQ in "${FASTQ_FILES[@]}"; do

    start_time=$(date +%s)

    # Extract sample name
    filename=$(basename "$FASTQ")
    sample="${filename}"
    for suffix in "_merged.fastq.gz" "_merged.fq.gz" "_merged.fastq" "_merged.fq"; do
        if [[ "$sample" == *"$suffix" ]]; then
            sample="${sample%"$suffix"}"
            break
        fi
    done

    echo "------------------------------------------------------------"
    echo "          Processing sample: $sample"
    echo "------------------------------------------------------------"


    # Identify reads containing the amplicon, trim both flanks, and convert to FASTA.
    echo "Selecting and extracting amplicon reads to FASTA..."

    seqkit amplicon --threads "$THREADS" -m 0 \
        -F "$FIVE_PRIME_FLANK" \
        -R "$REVERSE_AMPLICON" \
        -r "$(( ${#FIVE_PRIME_FLANK} + 1 )):-$(( ${#REVERSE_AMPLICON} + 1 ))" \
        "$FASTQ" | \
        seqkit fq2fa --threads "$THREADS" \
        > "${TRIMMED_DIR}/${sample}.fasta"


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
echo ""
echo ""
echo ""
