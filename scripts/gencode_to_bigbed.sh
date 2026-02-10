#!/usr/bin/env bash
set -euo pipefail

#############################################
# Convert GENCODE GTF → bigBed with gene symbols
# Author: Tomàs Montserrat / HERVarium project
#
# Usage:
#   ./gencode_to_bigbed.sh [gencode_release]
#
# This will:
#   - Download GENCODE v48 GTF + GRCh38 primary assembly FASTA
#   - Build chrom.sizes
#   - Convert GTF → genePred → BED
#   - Replace transcript IDs with gene symbols
#   - Build gencode.v48.genesymbols.bb (ready for IGV.js / Dash)
#
# Requirements: gtfToGenePred, genePredToBed, bedToBigBed, gawk, samtools, wget, gunzip
#############################################

# -------- Args --------
if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [gencode_release]" >&2
    exit 1
fi

REL="${1:-48}"
echo "[INFO] Using GENCODE release v${REL}"

# -------- Tool checks --------
for cmd in samtools sort cut gawk wget gunzip gtfToGenePred genePredToBed; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required tool '$cmd' not found in PATH." >&2
        exit 1
    fi
done

# Vendored UCSC binary (mandatory)
BEDTOBIGBED="./bin/bedToBigBed"
if [[ ! -x "$BEDTOBIGBED" ]]; then
    echo "ERROR: vendored bedToBigBed not found or not executable at $BEDTOBIGBED" >&2
    exit 1
fi

# -------- 1. Download inputs --------
GTFFILE="assets/gencode.v${REL}.primary_assembly.annotation.gtf"
FAFILE="assets/GRCh38.primary_assembly.genome.fa"

if [[ ! -f "$GTFFILE" ]]; then
    echo "[INFO] Downloading GENCODE GTF..."
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}/gencode.v${REL}.primary_assembly.annotation.gtf.gz"
    gunzip "gencode.v${REL}.primary_assembly.annotation.gtf.gz"
fi

if [[ ! -f "$FAFILE" ]]; then
    echo "[INFO] Downloading GRCh38 primary assembly FASTA..."
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}/GRCh38.primary_assembly.genome.fa.gz"
    gunzip "GRCh38.primary_assembly.genome.fa.gz"
fi

# -------- 2. Chromosome sizes --------
CHROMSIZES="${FAFILE}.chrom.sizes"
if [[ ! -f "$CHROMSIZES" ]]; then
    echo "[INFO] Computing chromosome sizes..."
    samtools faidx "$FAFILE"
    cut -f1,2 "${FAFILE}.fai" > "$CHROMSIZES"
fi

# -------- 3. GTF → genePred --------
GENEPRED="gencode.v${REL}.genepred"
if [[ ! -f "$GENEPRED" ]]; then
    echo "[INFO] Converting GTF → genePred..."
    gtfToGenePred -genePredExt "$GTFFILE" "$GENEPRED"
fi

# -------- 4. genePred → BED --------
BED="gencode.v${REL}.bed"
if [[ ! -f "$BED" ]]; then
    echo "[INFO] Converting genePred → BED..."
    genePredToBed "$GENEPRED" "$BED"
fi

# -------- 5. Replace transcript IDs with gene symbols --------
BED_SYMS="gencode.v${REL}.genesymbols.bed"
echo "[INFO] Replacing transcript IDs with gene symbols..."
gawk -F'\t' '
NR==FNR && $3=="transcript" {
    if (match($9, /transcript_id "([^"]+)"/, t) &&
        match($9, /gene_name "([^"]+)"/, g))
        map[t[1]] = g[1];
    next
}
NR>FNR {
    if ($4 in map) $4 = map[$4];
    print
}' "$GTFFILE" "$BED" > "$BED_SYMS"

# -------- 6. Sort BED --------
BED_SORTED="gencode.v${REL}.genesymbols.sorted.bed"
echo "[INFO] Sorting BED..."
sort -k1,1 -k2,2n "$BED_SYMS" > "$BED_SORTED"

# -------- 7. Build bigBed --------
BB_OUT="gencode.v${REL}.genesymbols.bb"
echo "[INFO] Building bigBed → $BB_OUT"
"$BEDTOBIGBED" "$BED_SORTED" "$CHROMSIZES" "$BB_OUT"

echo "[OK] Done."
echo "bigBed saved as: $BB_OUT"
