#!/usr/bin/env bash
# simplify_bed_names.sh
#
# Usage:
#   ./simplify_bed_names.sh HERV_domains_v1.bed HERV_domains_v1.simplified.bed
#
# It rewrites column 4 to: Subfamily|DomainType|Coverage(3dp)

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bed output.bed"
    exit 1
fi

infile="$1"
outfile="$2"

awk -F'\t' '
BEGIN{OFS="\t"}
# Pass-through track/header lines if any
/^track[[:space:]]/ || /^#/ {print; next}

{
  name = $4

  # Subfamily = everything before "-int" (if present)
  subfam = name
  if (match(name, /^(.*?)-int/, m)) {
      subfam = m[1]
  }

  # Split by pipes to get domain type (3rd) and coverage (6th)
  n = split(name, a, /\|/)
  domtype = (n >= 3 ? a[3] : "NA")
  covraw  = (n >= 6 ? a[6] : "")

  # Round coverage to 3 decimals if present
  cov = (covraw != "" ? sprintf("%.3f", covraw + 0) : "NA")

  # Replace column 4
  $4 = subfam "|" domtype "|" cov

  print
}' "$infile" > "$outfile"

echo "[OK] Simplified BED written to $outfile"


