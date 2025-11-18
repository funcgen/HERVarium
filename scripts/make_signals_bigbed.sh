#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_signals_bigbed.sh GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_signals.bed
#
# Output:
#   - assets/hg38.chrom.sizes
#   - assets/HERV_LTR_U3_R_U5_signals.clean.bed
#   - assets/HERV_LTR_U3_R_U5_signals.sorted.bed
#   - assets/HERV_LTR_U3_R_U5_signals.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <signals_bed>"
    echo "Example: $0 GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_signals.bed"
    exit 1
fi

GENOME_FA=$1
SIGNALS_BED=$2

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/HERV_LTR_U3_R_U5_signals.clean.bed"
SORTED_BED="assets/HERV_LTR_U3_R_U5_signals.sorted.bed"
BIGBED_OUT="assets/HERV_LTR_U3_R_U5_signals.bb"

# 0. Check required tools
for cmd in samtools ./bedToBigBed sort cut awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: '$cmd' not found in PATH. Please install it before running this script."
        exit 1
    fi
done

echo "1) Indexing FASTA with samtools faidx..."
samtools faidx "$GENOME_FA"

echo "   Extracting chromosome sizes -> $CHROM_SIZES"
cut -f1,2 "${GENOME_FA}.fai" > "$CHROM_SIZES"

echo "2) Cleaning signals BED (name -> TSS/Inr/..., integer score) -> $CLEAN_BED"
awk 'BEGIN{OFS="\t"}
{
    # Skip malformed / empty lines
    if (NF < 3) next

    # Ensure at least 6 columns (BED6)
    for (i = NF+1; i <= 6; i++) {
        if (i == 4) $i = "."
        else if (i == 5) $i = 0
        else if (i == 6) $i = "+"
    }

    # --- Rewrite name (col 4) ---
    # Original: seqname|SIG:TSS  or  seqname|SIG:Inr  etc.
    # We want:  TSS / Inr / DCE_SIII / PAS / CLEAVAGE ...
    split($4, a, "|")
    sigpart = (length(a[2]) ? a[2] : $4)   # SIG:TSS
    split(sigpart, b, ":")
    if (length(b[2])) {
        $4 = b[2]
    } else {
        $4 = sigpart
    }

    # --- Score (col 5) must be integer 0–1000 for bigBed ---
    score = $5
    if (score == "." || score == "" || score == "NA") {
        s = 0
    } else {
        s = int(score)   # your file already has small integers (1,2,...)
        if (s < 0)   s = 0
        if (s > 1000) s = 1000
    }
    $5 = s

    print
}' "$SIGNALS_BED" > "$CLEAN_BED"

echo "3) Sorting cleaned BED -> $SORTED_BED"
sort -k1,1 -k2,2n "$CLEAN_BED" > "$SORTED_BED"

echo "4) Converting to bigBed -> $BIGBED_OUT"
./bedToBigBed "$SORTED_BED" "$CHROM_SIZES" "$BIGBED_OUT"

echo "Done."
echo "Generated:"
echo "  - $CHROM_SIZES"
echo "  - $CLEAN_BED"
echo "  - $SORTED_BED"
echo "  - $BIGBED_OUT"
