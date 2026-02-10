#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_pbs_ppt_bigbed.sh GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_PBS_PPT.bed
#
# Output:
#   - assets/hg38.chrom.sizes
#   - assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.clean.bed
#   - assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.sorted.bed
#   - assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.bb

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome_fasta> <pbs_ppt_bed>" >&2
    echo "Example: $0 GRCh38.primary_assembly.genome.fa assets/HERV_LTR_U3_R_U5_PBS_PPT.bed" >&2
    exit 1
fi

GENOME_FA="$1"
PBS_PPT_BED="$2"

CHROM_SIZES="assets/hg38.chrom.sizes"
CLEAN_BED="assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.clean.bed"
SORTED_BED="assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.sorted.bed"
BIGBED_OUT="assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.bb"

# Ensure output dirs exist
mkdir -p "assets" "assets/ltr/segments"

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
        else if (i == 6) $i = "."   # safer default than "+"
    }

    # --- Rewrite name (col 4) to PBS / PPT only ---
    # Original: seqname|PBS or seqname|PPT
    split($4, a, "|")
    if (length(a[2])) {
        $4 = a[2]
    }

    # --- Score (col 5) must be integer 0–1000 for bigBed ---
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
}' "$PBS_PPT_BED" > "$CLEAN_BED"

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
