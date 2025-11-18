#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_pbs_ppt_bigbed.sh GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_PBS_PPT.bed
#
# Output:
#   - assets/hg38.chrom.sizes                      (reused by other scripts)
#   - assets/HERV_LTR_U3_R_U5_PBS_PPT.clean.bed
#   - assets/HERV_LTR_U3_R_U5_PBS_PPT.sorted.bed
#   - assets/HERV_LTR_U3_R_U5_PBS_PPT.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <pbs_ppt_bed>"
    echo "Example: $0 GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_PBS_PPT.bed"
    exit 1
fi

GENOME_FA=$1
PBS_PPT_BED=$2

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/HERV_LTR_U3_R_U5_PBS_PPT.clean.bed"
SORTED_BED="assets/HERV_LTR_U3_R_U5_PBS_PPT.sorted.bed"
BIGBED_OUT="assets/HERV_LTR_U3_R_U5_PBS_PPT.bb"

# 0. Check required tools
for cmd in samtools ./bedToBigBed sort cut awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: '$cmd' not found in PATH. Please install it before running this script."
        exit 1
    fi
done

echo "1) Indexing FASTA with samtools faidx (if needed)..."
samtools faidx "$GENOME_FA"

echo "   Extracting chromosome sizes -> $CHROM_SIZES"
cut -f1,2 "${GENOME_FA}.fai" > "$CHROM_SIZES"

echo "2) Cleaning PBS/PPT BED (name -> PBS/PPT, integer score) -> $CLEAN_BED"
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

    # --- Rewrite name (col 4) to PBS / PPT only ---
    # Original: seqname|PBS or seqname|PPT
    split($4, a, "|")
    if (length(a[2])) {
        $4 = a[2]   # PBS / PPT
    }

    # --- Score (col 5) must be integer 0–1000 for bigBed ---
    score = $5
    if (score == "." || score == "" || score == "NA") {
        s = 0
    } else {
        s = int(score)   # your file already has integers (10, 11, 12, ...)
        if (s < 0)   s = 0
        if (s > 1000) s = 1000
    }
    $5 = s

    print
}' "$PBS_PPT_BED" > "$CLEAN_BED"

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
