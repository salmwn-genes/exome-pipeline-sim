#!/usr/bin/env bash
# Exome-style variant calling pipeline (simulated data)
# FASTQ -> QC -> BWA-MEM -> sort -> MarkDuplicates -> HaplotypeCaller -> filter -> SnpEff
set -euo pipefail

REF=ref.fa
R1=sample_R1.fastq
R2=sample_R2.fastq
SAMPLE=sample1
OUT=results
mkdir -p "$OUT" qc

echo "[1/8] Simulating data"
[ -f "$REF" ] || python3 simulate.py

echo "[2/8] FastQC"
fastqc "$R1" "$R2" --outdir qc

echo "[3/8] Reference indexes"
bwa index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fa}.dict" || true

echo "[4/8] BWA-MEM alignment"
bwa mem -R "@RG\tID:sim1\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" "$REF" "$R1" "$R2" > "$OUT/aligned.sam"

echo "[5/8] Sort, mark duplicates, index"
samtools sort -o "$OUT/sorted.bam" "$OUT/aligned.sam"
gatk MarkDuplicates -I "$OUT/sorted.bam" -O "$OUT/dedup.bam" -M "$OUT/dup_metrics.txt"
samtools index "$OUT/dedup.bam"
samtools flagstat "$OUT/dedup.bam" > "$OUT/flagstat.txt"

echo "[6/8] Variant calling (GATK HaplotypeCaller)"
gatk HaplotypeCaller -R "$REF" -I "$OUT/dedup.bam" -O "$OUT/raw.vcf"

echo "[7/8] Filtering (QUAL>=30, DP>=10)"
bcftools filter -i 'QUAL>=30 && FORMAT/DP>=10' "$OUT/raw.vcf" -o "$OUT/filtered.vcf"

echo "[8/8] Annotation (SnpEff, synthetic gene model)"
snpEff ann -c snpEff.config chrSimDB "$OUT/filtered.vcf" > "$OUT/annotated.vcf"
bcftools query -f '%POS\t%REF\t%ALT\t%INFO/ANN\n' "$OUT/annotated.vcf" \
 | awk -F'\t' '{split($4,a,","); split(a[1],b,"|"); print $1,$2,$3,b[2],b[3],b[4],b[11]}' OFS='\t' \
 > "$OUT/annotation_table.tsv"

echo "[eval] Comparing with truth.txt"
bcftools query -f '%POS\t%REF\t%ALT\n' "$OUT/filtered.vcf" | sort > "$OUT/called.tsv"
awk '{print $2"\t"$3"\t"$4}' truth.txt | sort > "$OUT/truth_sorted.tsv"
echo "True positives:  $(comm -12 "$OUT/called.tsv" "$OUT/truth_sorted.tsv" | wc -l)"
echo "False positives: $(comm -23 "$OUT/called.tsv" "$OUT/truth_sorted.tsv" | wc -l)"
echo "False negatives: $(comm -13 "$OUT/called.tsv" "$OUT/truth_sorted.tsv" | wc -l)"
echo "Done. Results in $OUT/"
