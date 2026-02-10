#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_ltr_bigbed.sh assets/GRCh38.primary_assembly.genome.fa assets/ltr/ERVs_LTRs_merged_v4.simplified.bed
#
# Output:
#   - assets/hg38.chrom.sizes
#   - assets/ltr/ERV_ltr_merged.simplified.clean.bed
#   - assets/ltr/ERV_ltr_merged.simplified.sorted.bed
#   - assets/ltr/ERV_ltr_merged.simplified.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <ltr_bed>" >&2
    exit 1
fi

GENOME_FA="$1"
LTR_BED="$2"

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/ltr/ERV_ltr_merged.simplified.clean.bed"
SORTED_BED="assets/ltr/ERV_ltr_merged.simplified.sorted.bed"
BIGBED_OUT="assets/ltr/ERV_ltr_merged.simplified.bb"

# Ensure output dirs exist
mkdir -p "assets" "assets/ltr"

# -------- Tool checks --------
for cmd in samtools sort cut awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: '$cmd' not found in PATH." >&2
        exit 1
    fi
done

BEDTOBIGBED="./bin/bedToBigBed"
if [[ ! -x "$BEDTOBIGBED" ]]; then
    echo "ERROR: vendored bedToBigBed not found or not executable at $BEDTOBIGBED" >&2
    exit 1
fi

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
    for (i = NF+1; i <= 6; i++) {
        if (i == 4) $i = "."
        else if (i == 5) $i = 0
        else if (i == 6) $i = "."   # safer default than "+"
    }

    # Fix score column to integer 0–1000
    score = $5
    if (score == "." || score == "" || score == "NA") {
        s = 0
    } else {
        s = int(score)
        if (s < 0) s = 0
        if (s > 1000) s = 1000
    }
    $5 = s

    print
}' "$LTR_BED" > "$CLEAN_BED"

echo "3) Sorting cleaned BED -> $SORTED_BED"
sort -k1,1 -k2,2n "$CLEAN_BED" > "$SORTED_BED"

echo "4) Converting to bigBed -> $BIGBED_OUT"
"$BEDTOBIGBED" "$SORTED_BED" "$CHROM_SIZES" "$BIGBED_OUT"

echo "[OK] Done."
echo "Generated:"
echo "  - $CHROM_SIZES"
echo "  - $CLEAN_BED"
echo "  - $SORTED_BED"
echo "  - $BIGBED_OUT"
