#!/usr/bin/env bash
set -euo pipefail

# HERVarium: simplify BED "name" (col 4) to compact locus IDs
# - Internal:  PREFIX_pos_CHR_START_END_strand_S  ->  PREFIX_CHR_START_END_S
# - LTR:       PREFIX_merged_pos_CHR_START_END_strand_S -> PREFIX_CHR_START_END_S
#
# Usage:
#   ./simplify_internal_and_ltr_names.sh --type internal in.bed[.gz] out.bed[.gz]
#   ./simplify_internal_and_ltr_names.sh --type ltr      in.bed[.gz] out.bed[.gz]
#
# Notes:
# - CHR token supports non-"chr" contigs (e.g., KI270727.1, GL000220.1)
# - Only column 4 is modified; all other columns are preserved.
# - Prints a short summary (#lines processed / changed).

usage() {
  cat <<EOF
Usage: $0 --type {internal|ltr} <input.bed[.gz]> <output.bed[.gz]>

Examples:
  $0 --type internal assets/ERV_GyDB_v4_merged.bed assets/ERV_GyDB_v4_merged.simplified.bed
  $0 --type ltr      assets/ERV_ltr_v1_merged.bed   assets/ERV_ltr_v1_merged.simplified.bed
EOF
}

if [[ $# -lt 3 ]]; then
  usage; exit 1
fi

TYPE=""
IN=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --type) TYPE="${2:-}"; shift 2;;
    -h|--help) usage; exit 0;;
    -*)
      echo "Unknown option: $1" >&2; usage; exit 1;;
    *)
      if [[ -z "${IN:-}" ]]; then
        IN="$1"
      elif [[ -z "${OUT:-}" ]]; then
        OUT="$1"
      else
        echo "Too many positional arguments." >&2; usage; exit 1
      fi
      shift;;
  esac
done

if [[ -z "$TYPE" || -z "$IN" || -z "$OUT" ]]; then
  usage; exit 1
fi

if [[ "$TYPE" != "internal" && "$TYPE" != "ltr" ]]; then
  echo "ERROR: --type must be 'internal' or 'ltr'." >&2
  exit 1
fi

# Regex for gensub(): use [0-9]+ (portable) instead of \d+
if [[ "$TYPE" == "internal" ]]; then
  # ERVL-E-int_pos_chr1_41380_42285_strand_- -> ERVL-E-int_chr1_41380_42285_-
  REGEX='^(.+?)_pos_([^_]+)_([0-9]+)_([0-9]+)_strand_([+-])$'
  REPL='\\1_\\2_\\3_\\4_\\5'
else
  # MLT1K_merged_pos_chr1_21949_22344_strand_+ -> MLT1K_chr1_21949_22344_+
  REGEX='^(.+?)_merged_pos_([^_]+)_([0-9]+)_([0-9]+)_strand_([+-])$'
  REPL='\\1_\\2_\\3_\\4_\\5'
fi

# Readers/writers for gz or plain files
read_cmd()  { [[ "$1" =~ \.gz$ ]] && echo "gzip -cd --" || echo "cat --"; }
write_cmd() { [[ "$1" =~ \.gz$ ]] && echo "gzip -c >"  || echo "cat >"; }

READER=$(read_cmd "$IN")
WRITER=$(write_cmd "$OUT")

tmp_out="$(mktemp)"
trap 'rm -f "$tmp_out"' EXIT

# Use gawk (for gensub). If not present, try awk and hope it's GNU awk.
AWK_BIN="${AWK_BIN:-gawk}"

eval "$READER \"$IN\"" | "$AWK_BIN" -F'\t' -v OFS='\t' -v RGX="$REGEX" -v RPL="$REPL" '
BEGIN {changed=0; total=0}
{
  total++
  if (NF >= 4) {
    old=$4
    new=gensub(RGX, RPL, 1, old)
    if (new != old) changed++
    $4=new
  }
  print
}
END {
  printf("[INFO] Processed %d lines; changed %d names.\n", total, changed) > "/dev/stderr"
}' > "$tmp_out"

eval "$WRITER \"$OUT\"" < "$tmp_out"
echo "[OK] Wrote: $OUT"

