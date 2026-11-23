
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Backward-compatible positional form:
# merge_reads.sh INPUT_DIR [OUTPUT_DIR] [REPORT_DIR]
if [[ $# -gt 0 && "$1" != -* ]]; then
    INPUT_DIR="$1"
    OUTPUT_DIR="${2:-data/FASTQ_merged}"
    REPORT_DIR="${3:-results/fastp_reports}"
    exec python3 "$SCRIPT_DIR/merge_reads.py" \
        --input-dir "$INPUT_DIR" \
        --output-dir "$OUTPUT_DIR" \
        --report-dir "$REPORT_DIR"
fi

exec python3 "$SCRIPT_DIR/merge_reads.py" "$@"
