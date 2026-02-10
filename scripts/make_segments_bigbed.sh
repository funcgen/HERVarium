#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_segments_bigbed.sh GRCh38.primary_assembly.genome.fa HERV_LTR_U3_R_U5_segments_allconf.bed
#
# Output:
#   - assets/hg38.chrom.sizes
#   - assets/HERV_LTR_U3_R_U5_segments_allconf.clean.bed
#   - assets/HERV_LTR_U3_R_U5_segments_allconf.sorted.bed
#   - assets/HERV_LTR_U3_R_U5_segments_allconf.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <segments_bed>"
    echo "Example: $0 GRCh38.primary_assembly.genome.fa HERV_LTR_U3_R_U5_segments_allconf.bed"
    exit 1
fi

GENOME_FA=$1
SEGMENTS_BED=$2

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/ltr/segments/HERV_LTR_U3_R_U5_segments_allconf.clean.bed"
SORTED_BED="assets/ltr/segments/HERV_LTR_U3_R_U5_segments_allconf.sorted.bed"
BIGBED_OUT="assets/ltr/segments/HERV_LTR_U3_R_U5_segments_allconf.bb"

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

echo "2) Cleaning segments BED (name + score) -> $CLEAN_BED"
awk 'BEGIN{OFS="\t"}
{
    # Skip clearly malformed lines
    if (NF < 3) next

    # Ensure at least 6 columns (BED6)
    for (i = NF+1; i <= 6; i++) {
        if (i == 5) $i = 0
        else if (i == 6) $i = "+"
        else $i = "."
    }

    # --- Rewrite name (col 4) to U3 / U3|LOW_CONF / R / U5 ---

    # Original name looks like: seqname|U3|LOW_CONF
    # We keep token2 and optionally token3
    split($4, a, "|")
    newname = (length(a[2]) ? a[2] : $4)
    if (length(a[3])) {
        newname = newname "|" a[3]
    }
    $4 = newname

    # --- Fix score (col 5) to integer 0-1000 for bigBed ---

    score = $5
    if (score == "." || score == "" || score == "NA") {
        s = 0
    } else {
        s = int(score * 100)   # e.g. 1.05 -> 105, 4.83 -> 483
        if (s < 0)   s = 0
        if (s > 1000) s = 1000
    }
    $5 = s

    print
}' "$SEGMENTS_BED" > "$CLEAN_BED"

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

