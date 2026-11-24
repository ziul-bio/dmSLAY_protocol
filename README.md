# dmSLAY protocol

This repository provides a reproducible Pixi workflow for processing paired amplicon reads:
merge reads, trim amplicon flanks and adapters, count peptide sequences, and optionally
run differential-expression analysis with DESeq2.

## 1. Install the environment

Install [Pixi](https://pixi.sh/), then run these commands from the repository root:

```bash
pixi install
pixi run analyze --help
```

Pixi installs the tools and R/Python packages declared in `pixi.toml`. The committed
`pixi.lock` records the resolved versions for reproducible environments.

## 2. Prepare paired-read inputs

Place the FASTQ files in an input directory. Every sample must have both mates.
The workflow recognizes common naming patterns such as:

```text
sample_R1.fastq.gz   sample_R2.fastq.gz
sample_1.fastq.gz    sample_2.fastq.gz
sample_F.fastq.gz    sample_R.fastq.gz
sample_L001_R1_001.fastq.gz   sample_L001_R2_001.fastq.gz
```

If a mate is missing, merging stops and prints the sample and unmatched filename.

The 5-prime and 3-prime values are the flanking sequences immediately outside the
region of interest. They default to the example dataset values:

```text
5-prime: CTCCAGCTGCGGGTATCGGAGGAACC
3-prime: GTCGACCTGCAGGCATG
```

For amplicon selection, the 5-prime flank is passed to `seqkit amplicon -F`. The
3-prime flank is reverse-complemented for `seqkit amplicon -R`. By default, the
original flanks are also used as the left and right Flexbar adapter sequences.
Supply adapter FASTA files to override either default.

## 3. Run the full analysis

Run the complete workflow with:

```bash
pixi run analyze --input-dir /path/to/paired/reads
```

To override the flank defaults:

```bash
pixi run analyze --input-dir /path/to/paired/reads \
  --five-prime-sequence CTCCAGCTGCGGGTATCGGAGGAACC \
  --three-prime-sequence GTCGACCTGCAGGCATG
```

The workflow runs these stages in order:

- **Merge:** `fastp` merges each pair into `data/FASTQ_merged/` and writes HTML/JSON
  reports to `results/fastp_reports/`. Unmerged R1/R2 reads are also retained.
- **Trim:** `seqkit` selects the amplicon and `flexbar` trims adapters, writing FASTA
  files to `data/FASTA_trimmed/`.
- **Count:** translated peptide sequences are counted across samples and written to
  `results/counts/peptide_read_counts.txt`. By default, peptides with fewer than 10
  total reads are filtered out.
- **DESeq2 (optional):** runs only when you provide `--coldata-file` or `--groups`.
  Results are written to `results/deseq2/` by default.

The main workflow options are:

| Option | Purpose | Default |
| --- | --- | --- |
| `--input-dir DIR` | Paired FASTQ input directory | `data/FASTQ_raw/` |
| `--merged-dir DIR` | Merged FASTQ output | `data/FASTQ_merged/` |
| `--trimmed-dir DIR` | Trimmed FASTA output | `data/FASTA_trimmed/` |
| `--results-dir DIR` | Peptide-count output directory | `results/counts/` |
| `--report-dir DIR` | fastp HTML/JSON reports | `results/fastp_reports/` |
| `--threads N` | Worker threads for trimming tools | `32` |
| `--minimum-count N` | Minimum total reads per peptide | `10` |
| `--five-prime-sequence SEQ` | 5-prime flank | Example sequence above |
| `--three-prime-sequence SEQ` | 3-prime flank | Example sequence above |
| `--adapter-left FILE` | Left adapter FASTA override | 5-prime flank |
| `--adapter-right FILE` | Right adapter FASTA override | 3-prime flank |

Run `pixi run analyze --help` for the complete option list.

## 4. Configure optional DESeq2 analysis

DESeq2 requires a group for every sample. Provide groups inline or in a TSV/CSV file.

For inline groups, use comma- or semicolon-separated `SAMPLE=GROUP` entries:

```bash
pixi run analyze --input-dir /path/to/paired/reads \
  --groups "S1=Control,S2=Control,S3=Control,S4=Treated,S5=Treated,S6=Treated" \
  --deseq2-reference Control \
  --deseq2-treatment Treated
```

For a metadata file, the script recognizes these sample-column names:
`sample`, `sample_id`, `sample_name`, `Sample`, `SampleID`, or `name`. It recognizes
these group-column names: `group`, `condition`, `treatment`, or `phenotype`.

```tsv
sample	condition
S1	Control
S2	Control
S3	Control
S4	Treated
S5	Treated
S6	Treated
```

```bash
pixi run analyze --input-dir /path/to/paired/reads \
  --coldata-file /path/to/coldata.tsv \
  --deseq2-reference Control \
  --deseq2-treatment Treated
```

Every count-table sample must occur exactly once in the metadata. If no reference
group is specified, the first group is used. If no treatment is specified, each
non-reference group is compared against the reference. Set `--deseq2-output-dir`
to change the output location. The minimum count for DESeq2 is 30 by default.

## 5. Run individual stages

These tasks are useful when you want to run or repeat a single stage. Run them in
pipeline order, using the previous stage's output as the next stage's input.

```bash
pixi run merge --input-dir /path/to/paired/reads
pixi run trim -i data/FASTQ_merged \
  -5 CTCCAGCTGCGGGTATCGGAGGAACC \
  -3 GTCGACCTGCAGGCATG
pixi run count -i data/FASTA_trimmed -o results/counts
pixi run deseq2 \
  --count-file results/counts/peptide_read_counts.txt \
  --coldata-file data/coldata.tsv \
  --reference Control \
  --treatment Treated \
  --output-dir results/deseq2
```

## Troubleshooting installation

If Pixi is interrupted while installing or got stuck, rerun `pixi install`.

```bash
# First remove the generated workspace environment
rm -rf .pixi
rm -f pixi.lock

# And Then: On a network-backed cache, point Pixi to a local cache to reduce filesystem I/O:
export PIXI_CACHE_DIR="/tmp/pixi-cache-${USER}"
export PIXI_NO_PROGRESS=1

# and recreate it:
pixi install
```
