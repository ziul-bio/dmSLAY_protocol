#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

from pair_fastq import discover_pairs


def main():
    parser = argparse.ArgumentParser(description="Merge paired FASTQ reads with fastp")
    parser.add_argument("-i", "--input-dir", required=True, type=Path)
    parser.add_argument(
        "-o", "--output-dir", type=Path,
        help="Merged FASTQ output directory (default: results/<parent-name>/FASTQ_merged)",
    )
    parser.add_argument(
        "-r", "--report-dir", type=Path,
        help="fastp report directory (default: results/<parent-name>/fastp_reports)",
    )
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        parser.error(f"input directory does not exist: {args.input_dir}")

    output_root = Path(__file__).resolve().parent.parent / "results" / args.input_dir.parent.name
    args.output_dir = args.output_dir or output_root / "FASTQ_merged"
    args.report_dir = args.report_dir or output_root / "fastp_reports"

    pairs, unpaired = discover_pairs(args.input_dir)
    if unpaired:
        for sample, reads in unpaired:
            print(
                f"Sample {sample} has an unpaired FASTQ file: "
                f"{', '.join(path.name for path in reads)}. "
                "Find or add its matching mate."
            )
        raise SystemExit(1)
    if not pairs:
        parser.error("no paired FASTQ files found")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    print("")
    print("************************************************************")
    print("          Merging paired FASTQ reads with fastp"             )
    print("************************************************************")

    for sample, read1, read2 in pairs:
        print("")
        print("-" * 60)
        print(f"              Merging sample: {sample}")
        print("-" * 60)
        print(f"R1: {read1}")
        print(f"R2: {read2}")

        command = [
            "fastp",
            "--in1", str(read1),
            "--in2", str(read2),
            "--merge",
            "--out1", str(args.output_dir / f"{sample}_R1_unmerged.fastq.gz"),
            "--out2", str(args.output_dir / f"{sample}_R2_unmerged.fastq.gz"),
            "--merged_out", str(args.output_dir / f"{sample}_merged.fastq.gz"),
            "--json", str(args.report_dir / f"{sample}.json"),
            "--html", str(args.report_dir / f"{sample}.html"),
        ]
        subprocess.run(command, check=True)

    print("-" * 60)
    print(f"Processing complete. Total samples processed: {len(pairs)}")
    print("\n" * 3, end="")  # adds three blank lines


if __name__ == "__main__":
    main()
