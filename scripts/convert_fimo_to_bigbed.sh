#!/usr/bin/env bash
# convert_bed_to_bigbed.sh
#
# Usage:
#   ./convert_bed_to_bigbed.sh input.bed genome.fa.fai
#
# Output:
#   - input.sorted.bed
#   - input.fixed.bed (score + name simplified)
#   - input.bb (bigBed)
#   - hg38.chrom.sizes

set -euo pipefail

# -------- Input check --------
if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 <input.bed> <genome.fai>" >&2
  exit 1
fi

BED_IN="$1"
FAI="$2"

BED_SORTED="${BED_IN%.bed}.sorted.bed"
BED_FIXED="${BED_IN%.bed}.fixed.bed"
CHROM_SIZES="hg38.chrom.sizes"
BB_OUT="${BED_IN%.bed}.bb"

# -------- Tool checks --------
for cmd in sort cut awk; do
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

# -------- 1. Sort BED --------
echo "[INFO] Sorting BED file..."
sort -k1,1 -k2,2n "$BED_IN" > "$BED_SORTED"

# -------- 2. Fix column 4 (simplify name) and column 5 (score as int 0–1000) --------
echo "[INFO] Simplifying column 4 and fixing column 5..."
awk 'BEGIN{OFS="\t"}
{
  # Example name:
  # MER34_merged_pos_GL000008.2_10110_10225_strand_+_MA2331.1_ZNF157
  name = $4
  split(name, parts, "_")
  subfamily = parts[1]

  # Extract TF from last token
  if (match(name, /_([^_]+)$/, m)) {
    tf = m[1]
  } else {
    tf = "NA"
  }

  $4 = subfamily "|" tf

  # Fix score → integer 0–1000
  score = int($5 + 0.5)
  if (score < 0) score = 0
  if (score > 1000) score = 1000
  $5 = score

  print
}' "$BED_SORTED" > "$BED_FIXED"

# -------- 3. Make chrom.sizes from .fai --------
echo "[INFO] Generating chrom.sizes from FASTA index..."
cut -f1,2 "$FAI" > "$CHROM_SIZES"

# -------- 4. Convert to bigBed --------
echo "[INFO] Converting to bigBed..."
"$BEDTOBIGBED" "$BED_FIXED" "$CHROM_SIZES" "$BB_OUT"

echo "[OK] Done. bigBed saved as: $BB_OUT"
