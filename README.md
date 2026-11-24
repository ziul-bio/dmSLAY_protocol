# dmSLAY protocol

This repository provides a reproducible Pixi workflow for merging paired amplicon reads, trimming their amplicon flanks/adapters, and counting peptide sequences.

## Requirements

Install [Pixi](https://pixi.sh/) and run commands from the repository root:

```bash
pixi install
pixi run analyze --input-dir /path/to/paired/reads \
  --five-prime-sequence CTCCAGCTGCGGGTATCGGAGGAACC \
  --three-prime-sequence GTCGACCTGCAGGCATG
```

The 5-prime and 3-prime values are the flanking sequences immediately outside the region of interest. They default to the example dataset values:

```text
5-prime: CTCCAGCTGCGGGTATCGGAGGAACC
3-prime: GTCGACCTGCAGGCATG
```

For amplicon selection, the 5-prime flank is used as `seqkit amplicon -F`, while the 3-prime flank is reverse-complemented and used as `-R`. For right-side trimming, Flexbar receives the original 3-prime flank.

## Input

All reads must be paired. Files may use common naming patterns, including:

- `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
- `sample_1.fastq.gz` and `sample_2.fastq.gz`
- `sample_F.fastq.gz` and `sample_R.fastq.gz`
- `sample_L001_R1_001.fastq.gz` and `sample_L001_R2_001.fastq.gz`

Missing mates and duplicate read roles are reported as errors rather than silently skipped.

## Options

```text
--input-dir DIR          Directory containing paired reads
--merged-dir DIR         Merged FASTQ output
--trimmed-dir DIR        Trimmed FASTA output
--results-dir DIR        Count table output
--report-dir DIR         fastp HTML/JSON reports
--adapter-left FILE      Optional left adapter FASTA override
--adapter-right FILE     Optional right adapter FASTA override
--five-prime-sequence SEQ   5-prime/forward amplicon flank
--three-prime-sequence SEQ  3-prime/reverse amplicon flank
--threads N              Sequencing worker threads
--minimum-count N        Minimum peptide read count (default: 10)
--coldata-file FILE      TSV/CSV metadata with sample and group/condition columns
--groups SPEC            Inline groups: SAMPLE=GROUP,SAMPLE=GROUP
--deseq2-reference GROUP DESeq2 reference group
--deseq2-treatment GROUP DESeq2 treatment group(s), comma-separated
--deseq2-output-dir DIR  DESeq2 output directory
--deseq2-minimum-count N Minimum peptide count for DESeq2 (default: 30)
```

Adapter trimming is skipped when adapter files are not supplied or do not exist.

## DESeq2 grouping

DESeq2 is run after the count table is generated when `--coldata-file` or `--groups`
is supplied.

The simplest option is an inline group specification:

```bash
pixi run analyze --input-dir /path/to/paired/reads \
  --groups "S1=Control,S2=Control,S3=Control,S4=Treated,S5=Treated,S6=Treated" \
  --deseq2-reference Control \
  --deseq2-treatment Treated
```

For larger experiments, provide a TSV or CSV file. Recognized sample columns are
`sample`, `sample_id`, `sample_name`, `Sample`, `SampleID`, or `name`; recognized
group columns are `group`, `condition`, `treatment`, or `phenotype`.

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

Every count-table sample must occur exactly once in the metadata. If no treatment
is supplied, DESeq2 compares every non-reference group against the reference.
Results and volcano plots are written to `results/deseq2/` by default.

## Direct tasks

```bash
pixi run merge --input-dir path/to/FASTQ
pixi run trim --input-dir data/FASTQ_merged \
  -5 CTCCAGCTGCGGGTATCGGAGGAACC \
  -3 GTCGACCTGCAGGCATG
pixi run count --input-dir data/FASTA_trimmed --output-dir results/counts

pixi run deseq2 \
  --count-file results/counts/peptide_read_counts.txt \
  --coldata-file data/coldata.tsv \
  --reference Control \
  --treatment Treated \
  --output-dir results/deseq2
```
