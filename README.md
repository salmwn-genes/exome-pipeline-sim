# Exome-style variant calling pipeline (simulated data)

An end-to-end germline variant calling workflow, from FASTQ to annotated VCF,
built to practice and document the standard GATK-based approach.

> **Note:** All data in this repository is **simulated**. The reference is a
> random 50 kb sequence (`chrSim`), the reads are synthetic paired-end reads
> (~100x), and the gene model is synthetic. Results demonstrate that the
> pipeline works; they say nothing about performance on real human exomes.

## Workflow

FASTQ → FastQC → BWA-MEM → samtools sort → GATK MarkDuplicates →
GATK HaplotypeCaller → bcftools filter → SnpEff annotation

## Tools

| Step | Tool | Version |
|---|---|---|
| QC | FastQC | 0.13.0 |
| Alignment | BWA-MEM | 0.7.19 |
| BAM processing | SAMtools | 1.25 |
| Duplicates, calling | GATK | 4.6.2.0 |
| Filtering | bcftools | 1.25 |
| Annotation | SnpEff | (conda, Java 21) |
| Visualization | IGV Desktop | |

Environment: Ubuntu on WSL2, Conda (Miniforge).

## Data simulation

`simulate.py` generates:
- a random 50 kb reference (`ref.fa`)
- 20 heterozygous SNPs (`truth.txt`)
- 25,000 paired-end read pairs (2x100 bp, 0.5% random base error)

`make_gff.py` defines 6 synthetic single-exon genes for annotation.

## Run

```bash
conda activate exome
./run_pipeline.sh
```

## Results

**Alignment (flagstat):** 50,000 reads, 100% mapped, 100% properly paired,
72 duplicates (0.14%).

**Variant calling vs. truth:**

| Metric | Value |
|---|---|
| True positives | 20 / 20 |
| False positives | 0 |
| False negatives | 0 |

All calls were heterozygous (0/1) with allele depths near 50/50.

**Annotation (SnpEff, synthetic genes):** 9 missense, 6 synonymous,
4 upstream/downstream, 1 intergenic. See `results/annotation_table.tsv`.

**IGV check:** variant at chrSim:1196 shows ~half of reads carrying the
alternate allele, while scattered single-read mismatches are sequencing
errors correctly ignored by the caller (`igv_chrSim_1196.png`).

## Limitations

- Simulated data: no repeats, no systematic errors, no real gene structure.
- Perfect accuracy is expected here and would not be on real exomes.
- No BQSR, no VQSR/hard filtering beyond QUAL/DP, no population database
  filtering (gnomAD) and no clinical interpretation (ACMG).

## Next steps

- Run on a real public exome subset against GRCh38
- Annotate with a human database (SnpEff/VEP), filter with gnomAD, check ClinVar
- Practice ACMG-based variant interpretation
