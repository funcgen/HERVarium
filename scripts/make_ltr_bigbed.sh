#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_ltr_bigbed.sh GRCh38.primary_assembly.genome.fa assets/ERV_ltr_v2_merged.simplified.bed
#
# Output:
#   - assets/hg38.chrom.sizes
#   - assets/ERV_ltr_v2_merged.simplified.clean.bed
#   - assets/ERV_ltr_v2_merged.simplified.sorted.bed
#   - assets/ERV_ltr_v2_merged.simplified.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <ltr_bed>"
    echo "Example: $0 GRCh38.primary_assembly.genome.fa assets/ERV_ltr_v2_merged.simplified.bed"
    exit 1
fi

GENOME_FA=$1
LTR_BED=$2

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/ltr/ERV_ltr_merged.simplified.clean.bed"
SORTED_BED="assets/ltr/ERV_ltr_merged.simplified.sorted.bed"
BIGBED_OUT="assets/ltr/ERV_ltr_merged.simplified.bb"

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

echo "2) Cleaning LTR BED (ensure BED6 + integer score) -> $CLEAN_BED"
awk 'BEGIN{OFS="\t"}
{
    # Skip malformed or empty lines
    if (NF < 3) next

    # Ensure at least 6 columns (BED6)
    # col1: chrom, col2: start, col3: end
    # col4: name (keep as-is or "." if missing)
    # col5: score (int 0–1000 for bigBed)
    # col6: strand ("+" / "-" / "." if missing)
    for (i = NF+1; i <= 6; i++) {
        if (i == 4) $i = "."
        else if (i == 5) $i = 0
        else if (i == 6) $i = "+"
    }

    # Fix score column to integer 0–1000
    score = $5
    if (score == "." || score == "" || score == "NA") {
        s = 0
    } else {
        s = int(score)        # if it was already integer; if float, this truncates
        if (s < 0)   s = 0
        if (s > 1000) s = 1000
    }
    $5 = s

    print
}' "$LTR_BED" > "$CLEAN_BED"

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

