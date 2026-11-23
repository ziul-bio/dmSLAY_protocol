#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


READ_EXTENSIONS = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)


def strip_read_extension(name: str) -> str:
    lower = name.lower()
    for suffix in READ_EXTENSIONS:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def split_paired_name(path: Path):
    stem = strip_read_extension(path.name)

    named = re.match(r"^(.*?)[._-](R1|R2|F|R)(?:[._-]?[0-9]+)?$", stem, re.IGNORECASE)
    if named:
        sample, role = named.groups()
        role = role.upper()
        return sample, ("1" if role in ("R1", "F") else "2")

    numeric = re.match(r"^(.*?)[._-]([12])(?:[._-]?[0-9]+)?$", stem)
    if numeric:
        sample, role = numeric.groups()
        return sample, role

    trailing_role = re.match(r"^(.*?)(?:[._-]?)(F|R)$", stem, re.IGNORECASE)
    if trailing_role:
        sample, role = trailing_role.groups()
        sample = sample.rstrip("._-")
        return sample, ("1" if role.upper() == "F" else "2")

    return None


def discover_pairs(input_dir: Path):
    pairs = {}
    unpaired = []

    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        split = split_paired_name(path)
        if split is None:
            continue
        sample, role = split
        pairs.setdefault(sample, {})[role] = path

    complete = []
    for sample in sorted(pairs):
        reads = pairs[sample]
        if "1" in reads and "2" in reads:
            complete.append((sample, reads["1"], reads["2"]))
        else:
            unpaired.append((sample, [reads[role] for role in sorted(reads)]))

    return complete, unpaired


def main():
    parser = argparse.ArgumentParser(
        description="Discover paired FASTQ files using common Illumina-style naming patterns."
    )
    parser.add_argument("input_dir", type=Path)
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        parser.error(f"input directory does not exist: {args.input_dir}")

    pairs, unpaired = discover_pairs(args.input_dir)
    if not pairs:
        parser.error("no paired FASTQ files found")
    if unpaired:
        parser.error(
            "missing mate for: "
            + ", ".join(
                f"{sample} ({', '.join(str(path.name) for path in reads)})"
                for sample, reads in unpaired
            )
        )

    for sample, read1, read2 in pairs:
        print(f"{sample}\t{read1}\t{read2}")


if __name__ == "__main__":
    main()
