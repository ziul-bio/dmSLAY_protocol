#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


FIELDS = (
    ("before_filtering", "total_reads"),
    ("before_filtering", "total_bases"),
    ("before_filtering", "q20_rate"),
    ("before_filtering", "q30_rate"),
    ("after_filtering", "total_reads"),
    ("after_filtering", "total_bases"),
    ("after_filtering", "q20_rate"),
    ("after_filtering", "q30_rate"),
    ("filtering_result", "passed_filter_reads"),
    ("filtering_result", "corrected_reads"),
    ("filtering_result", "low_quality_reads"),
    ("filtering_result", "too_many_N_reads"),
    ("filtering_result", "too_short_reads"),
    ("filtering_result", "too_long_reads"),
    ("duplication", "rate"),
    ("insert_size", "peak"),
    ("insert_size", "unknown"),
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract summary statistics from fastp JSON reports."
    )
    parser.add_argument("--input_dir", "-i", type=Path, required=True)
    parser.add_argument("--output", "-o", type=Path, required=True)
    args = parser.parse_args()

    reports = sorted(args.input_dir.glob("S*.json"))
    if not reports:
        raise FileNotFoundError(f"No fastp JSON reports found in {args.input_dir}")

    headers = ["sample"] + [f"{section}_{field}" for section, field in FIELDS]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        writer.writerow(headers)
        for report in reports:
            with report.open() as input_file:
                data = json.load(input_file)
            row = [report.stem]
            for section, field in FIELDS:
                try:
                    row.append(data[section][field])
                except KeyError as error:
                    raise KeyError(
                        f"Missing fastp field {section}.{field} in {report}"
                    ) from error
            writer.writerow(row)


if __name__ == "__main__":
    main()
