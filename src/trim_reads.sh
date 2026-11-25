
#!/usr/bin/env bash

set -euo pipefail

# Usage:
# bash scripts/trim_reads.sh \
#     -i data/FASTQ_merged \
#     -o data/FASTA_trimmed \
#     -t 32 \
#     -5 CTCCAGCTGCGGGTATCGGAGGAACC \
#     -3 GTCGACCTGCAGGCATG


# Default parameters
INPUT_DIR=""
TRIMMED_DIR="data/FASTA_trimmed"

ADAPTER_LEFT=""
ADAPTER_RIGHT=""

THREADS=50

FIVE_PRIME_FLANK="CTCCAGCTGCGGGTATCGGAGGAACC"
THREE_PRIME_FLANK="GTCGACCTGCAGGCATG"


# Parse arguments
while getopts ":i:o:a:b:t:f:r:5:3:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) TRIMMED_DIR="$OPTARG" ;;
        a) ADAPTER_LEFT="$OPTARG" ;;
        b) ADAPTER_RIGHT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f|5) FIVE_PRIME_FLANK="$OPTARG" ;;
        r|3) THREE_PRIME_FLANK="$OPTARG" ;;
        h)
            echo "Usage: $0 -i INPUT_DIR [options]"
            echo ""
            echo "Options:"
            echo "  -i  Input FASTQ directory (required)"
            echo "  -o  Output FASTA directory (default: $TRIMMED_DIR)"
            echo "  -a  Optional left adapter file (overrides -5)"
            echo "  -b  Optional right adapter file (overrides -3)"
            echo "  -t  Number of threads (default: $THREADS)"
            echo "  -5  5-prime (forward) flanking sequence"
            echo "  -3  3-prime (reverse) flanking sequence"
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


# Create output directory
mkdir -p "$TRIMMED_DIR"


# Create temporary directory and configure cleanup
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT



normalize_sequence() { # Normalize a DNA sequence to uppercase and remove whitespace
    printf '%s' "$1" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]'
}



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



FIVE_PRIME_FLANK="$(normalize_sequence "$FIVE_PRIME_FLANK")"
THREE_PRIME_FLANK="$(normalize_sequence "$THREE_PRIME_FLANK")"
REVERSE_AMPLICON="$(reverse_complement "$THREE_PRIME_FLANK")"



# Flexbar trims in read orientation. Use the 5' flank on the left and the
# original 3' flank on the right; only the seqkit -R sequence is reversed.
printf '>five_prime_flank\n%s\n' "$FIVE_PRIME_FLANK" > "${TEMP_DIR}/five_prime_flank.fasta"
printf '>three_prime_flank\n%s\n' "$THREE_PRIME_FLANK" > "${TEMP_DIR}/three_prime_flank.fasta"

if [[ -z "$ADAPTER_LEFT" ]]; then
    ADAPTER_LEFT="${TEMP_DIR}/five_prime_flank.fasta"
fi
if [[ -z "$ADAPTER_RIGHT" ]]; then
    ADAPTER_RIGHT="${TEMP_DIR}/three_prime_flank.fasta"
fi

if [[ -f "$ADAPTER_LEFT" ]]; then
    HAS_ADAPTER_LEFT=true
else
    HAS_ADAPTER_LEFT=false
fi
if [[ -f "$ADAPTER_RIGHT" ]]; then
    HAS_ADAPTER_RIGHT=true
else
    HAS_ADAPTER_RIGHT=false
fi


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
echo "                    Read Trimming"
echo "************************************************************"
echo "Input directory: $INPUT_DIR"
echo "Temporary directory: $TEMP_DIR"
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


    # Step 1: Filter reads by amplicon and convert to FASTA
    echo "Filtering reads by amplicon and converting to FASTA..."

    seqkit amplicon --threads "$THREADS" -m 2 \
        -F "$FIVE_PRIME_FLANK" \
        -R "$REVERSE_AMPLICON" \
        "$FASTQ" | \
        seqkit fq2fa --threads "$THREADS" \
        > "${TEMP_DIR}/${sample}.fasta"


    CURRENT="${TEMP_DIR}/${sample}.fasta"

    if [[ "$HAS_ADAPTER_LEFT" == true ]]; then
        echo "Trimming left adapter..."
        flexbar --threads "$THREADS" -ae 0.1 \
            --reads "$CURRENT" \
            --adapters "$ADAPTER_LEFT" \
            --adapter-trim-end LEFT \
            --target "${TEMP_DIR}/${sample}_left_tmp"
        CURRENT="${TEMP_DIR}/${sample}_left_tmp.fasta"
        if [[ ! -f "$CURRENT" ]]; then
            echo "Error: Flexbar did not generate the left-trimmed FASTA for $sample."
            exit 1
        fi
    fi

    if [[ "$HAS_ADAPTER_RIGHT" == true ]]; then
        echo "Trimming right adapter..."
        flexbar --threads "$THREADS" -ae 0.1 \
            --reads "$CURRENT" \
            --adapters "$ADAPTER_RIGHT" \
            --adapter-trim-end RIGHT \
            --target "${TRIMMED_DIR}/${sample}"
    elif [[ "$CURRENT" != "${TRIMMED_DIR}/${sample}.fasta" ]]; then
        cp "$CURRENT" "${TRIMMED_DIR}/${sample}.fasta"
    fi


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
