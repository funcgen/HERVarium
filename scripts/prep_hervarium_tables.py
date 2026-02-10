#!/usr/bin/env python3
import re, json
from pathlib import Path
import pandas as pd
import numpy as np

ASSETS = Path("assets")
OUTDIR = ASSETS / "precomputed"
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================
# Helpers
# =========================

# ERV locus normalizer:
#  "<prefix>_pos_<chrom>_<start>_<end>_strand_<+|->"
#   → "<prefix>_<chrom>_<start>_<end>_<+|->"
_LOCUS_NORM_RE = re.compile(
    r'^(?P<prefix>.+?)_pos_(?P<chrom>[^_]+)_(?P<start>\d+)_(?P<end>\d+)_strand_(?P<strand>[+-])$'
)
def normalize_internal_name(s: str):
    """Normalize coordinate-style names; preserve NA; leave enumerated IDs intact."""
    if pd.isna(s):
        return pd.NA
    s = str(s)
    m = _LOCUS_NORM_RE.match(s)
    return s if not m else f"{m.group('prefix')}_{m.group('chrom')}_{m.group('start')}_{m.group('end')}_{m.group('strand')}"

def strip_int_suffix(x: str) -> str:
    """Remove trailing '-int' or '_int' (case-insensitive)."""
    if pd.isna(x):
        return x
    s = str(x).strip()
    return re.sub(r'([_-])int$', '', s, flags=re.IGNORECASE)

_CHROM_SPAN_RE = re.compile(r'_(?P<chrom>[^_]+)_(?P<start>\d+)_(?P<end>\d+)(?:_[+-])?$')
def locus_to_interval(locus: str):
    """Extract (chrom, start, end) from a normalized locus name or enumerated (no-op)."""
    if pd.isna(locus):
        return (np.nan, np.nan, np.nan)
    m = _CHROM_SPAN_RE.search(str(locus))
    if not m:
        return (np.nan, np.nan, np.nan)
    chrom, start, end = m.group('chrom'), int(m.group('start')), int(m.group('end'))
    if end < start:
        start, end = end, start
    return (chrom, start, end)

# =========================
# Load inputs
# =========================

ltr_tsv = ASSETS / "ltr/LTR_fully_annotated.tsv"
int_full_tsv = ASSETS / "internals/INTERNAL_fully_annotated.tsv"
dom_tsv = ASSETS / "internals/HERV_loci_annotated_domains_v6.tsv"

df_ltr_raw = pd.read_csv(ltr_tsv, sep="\t", low_memory=False)
df_int_full = pd.read_csv(int_full_tsv, sep="\t", low_memory=False)
df_dom_raw = pd.read_csv(dom_tsv, sep="\t", low_memory=False)

# =========================
# Build coord → enumerated map from INTERNAL (source of truth)
# =========================
# Case 1 (old): INTERNAL has both 'internal_name' (enumerated ID) and 'name' (coord-style)
# Case 2 (new): INTERNAL only has 'name' → we use it as both ID and coord-style
if "internal_name" in df_int_full.columns:
    # Old style: enumerated + coord-based
    int_map = df_int_full[["internal_name", "name"]].copy()
    int_map.rename(columns={"internal_name": "enum_id"}, inplace=True)
else:
    # New style: only coord-based name; treat it as canonical ID as well
    if "name" not in df_int_full.columns:
        raise ValueError("INTERNAL_fully_annotated.tsv must contain a 'name' column")
    int_map = df_int_full[["name"]].copy()
    int_map.rename(columns={"name": "enum_id"}, inplace=True)

# Build normalized locus key from the coord-style name (which is *always* 'name' in the TSV)
# For the new format, enum_id == name == coord key; for the old format, coord key is 'name'
int_map["Locus_norm"] = int_map["enum_id"].astype(str).map(normalize_internal_name)

int_map = (
    int_map.dropna(subset=["Locus_norm"])
           .drop_duplicates(subset=["Locus_norm"])
)




# =========================
# LTR table
# =========================
df_ltr = df_ltr_raw.copy()

keep_ltr = [
    "sequence_name", "chrom", "start", "end", "strand",
    "RepeatName", "ltr_ClassFamily", "LTR_length", "type",
    "internal_name", "n_motifs", "capped_motifs",
    "load_category", "dist_to_tss", "nearest_gene", "tss_bin"
]
df_ltr = df_ltr[[c for c in keep_ltr if c in df_ltr.columns]].copy()

# Normalize sequence_name only (visual column)
pat = r'^(.+?)_merged_pos_([^_]+)_(\d+)_(\d+)_strand_([+-])$'
if "sequence_name" in df_ltr.columns:
    df_ltr["sequence_name"] = (
        df_ltr["sequence_name"].astype(str)
        .str.replace(pat, r"\1_\2_\3_\4_\5", regex=True)
        .str.replace("_merged_pos_", "_", regex=False)
        .str.replace("_strand_", "_", regex=False)
    )

df_ltr.rename(columns={
    "RepeatName": "subfamily",
    "start": "ltr_start", "end": "ltr_end", "strand": "ltr_strand",
    "LTR_length": "ltr_len",
}, inplace=True)

# Keep enumerated 'internal_name' as-is; derive int_subfamily for filtering
if "internal_name" in df_ltr.columns:
    df_ltr["int_subfamily"] = (
        df_ltr["internal_name"]
        .astype("string")
        .str.extract(r'^(.*?)-int', expand=False)
        .map(lambda s: strip_int_suffix(s) if pd.notna(s) else pd.NA)
    )

# Types
for c in ["subfamily", "ltr_ClassFamily", "ltr_strand", "type",
          "int_subfamily", "nearest_gene", "tss_bin", "internal_name"]:
    if c in df_ltr.columns:
        df_ltr[c] = df_ltr[c].astype("category")

# Save LTR parquet + meta
df_ltr.to_parquet(OUTDIR / "ltr.parquet", index=False)
ltr_meta = {
    "min_len":  int(np.nanmin(df_ltr["ltr_len"])) if "ltr_len" in df_ltr else 0,
    "max_len":  int(np.nanmax(df_ltr["ltr_len"])) if "ltr_len" in df_ltr else 0,
    "min_dist": int(np.nanmin(df_ltr["dist_to_tss"])) if "dist_to_tss" in df_ltr else 0,
    "max_dist": int(np.nanmax(df_ltr["dist_to_tss"])) if "dist_to_tss" in df_ltr else 0,
    "ltr_subfamilies": sorted(df_ltr["subfamily"].dropna().astype(str).unique().tolist()) if "subfamily" in df_ltr else [],
    "ltr_types":       sorted(df_ltr["type"].dropna().astype(str).unique().tolist()) if "type" in df_ltr else [],
    "ltr_int_subf":    sorted(df_ltr["int_subfamily"].dropna().astype(str).unique().tolist()) if "int_subfamily" in df_ltr else [],
}
(OUTDIR / "ltr_meta.json").write_text(json.dumps(ltr_meta, indent=2))


# =========================
# U3/R/U5 + PBS/PPT + promoter/PAS signals table
# =========================

def _load_u3r_bed(path: Path, source: str) -> pd.DataFrame:
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        dtype={"chrom": "string"},
        low_memory=False,
    )
    df["source"] = source
    return df

u3r_pbs_ppt_bed = ASSETS / "ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.bed"
u3r_segments_bed = ASSETS / "ltr/segments/HERV_LTR_U3_R_U5_segments_allconf.bed"
u3r_signals_bed = ASSETS / "ltr/segments/HERV_LTR_U3_R_U5_signals.bed"

df_pbs_ppt = _load_u3r_bed(u3r_pbs_ppt_bed, "PBS_PPT")
df_segments = _load_u3r_bed(u3r_segments_bed, "SEGMENT")
df_signals = _load_u3r_bed(u3r_signals_bed, "SIGNAL")

df_u3r_all = pd.concat([df_segments, df_pbs_ppt, df_signals], ignore_index=True)

# Parse the "name" field: "<sequence_name_raw>|TOKEN1[|TOKEN2]"
parts = df_u3r_all["name"].astype("string").str.split("|", n=2, expand=True)
df_u3r_all["sequence_name_raw"] = parts[0]
df_u3r_all["token1"] = parts[1]
df_u3r_all["token2"] = parts[2] if parts.shape[1] > 2 else pd.NA

def _classify_u3r(row):
    src = row["source"]
    t1 = row["token1"]
    t2 = row["token2"]

    if src == "SEGMENT":
        feature_class = "SEGMENT"
        feature = t1
        conf = t2 if pd.notna(t2) and str(t2) != "" else "HIGH_CONF"
    elif src == "PBS_PPT":
        feature_class = "PBS_PPT"
        feature = t1
        conf = pd.NA
    else:  # SIGNAL
        feature_class = "SIGNAL"
        if isinstance(t1, str) and t1.startswith("SIG:"):
            feature = t1.split(":", 1)[1]
        else:
            feature = t1
        conf = pd.NA

    return pd.Series(
        {
            "feature_class": feature_class,
            "feature": feature,
            "segment_conf": conf,
        }
    )

df_u3r_all = df_u3r_all.join(df_u3r_all.apply(_classify_u3r, axis=1))

# Normalise sequence_name so it matches df_ltr["sequence_name"]
df_u3r_all["sequence_name"] = (
    df_u3r_all["sequence_name_raw"]
    .astype("string")
    .str.replace(pat, r"\1_\2_\3_\4_\5", regex=True)
    .str.replace("_merged_pos_", "_", regex=False)
    .str.replace("_strand_", "_", regex=False)
)

# Attach LTR-level metadata from df_ltr
ltr_merge_cols = [
    c
    for c in [
        "sequence_name",
        "subfamily",
        "ltr_ClassFamily",
        "type",
        "internal_name",
        "int_subfamily",
        "ltr_start",
        "ltr_end",
        "ltr_len",
        "ltr_strand",
        "dist_to_tss",
        "nearest_gene",
        "tss_bin",
    ]
    if c in df_ltr.columns
]
df_ltr_small = df_ltr[ltr_merge_cols].copy()

df_u3r = df_u3r_all.merge(df_ltr_small, on="sequence_name", how="left")

# Feature length and pretty locus for display
df_u3r["feature_len"] = df_u3r["end"] - df_u3r["start"]
df_u3r["feature_locus"] = (
    df_u3r["chrom"].astype("string")
    + ":"
    + df_u3r["start"].astype("Int64").astype("string")
    + "-"
    + df_u3r["end"].astype("Int64").astype("string")
    + " ("
    + df_u3r["strand"].astype("string")
    + ")"
)

# Tidy dtypes
for c in [
    "feature_class",
    "feature",
    "segment_conf",
    "subfamily",
    "ltr_ClassFamily",
    "type",
    "int_subfamily",
    "nearest_gene",
    "tss_bin",
    "sequence_name",
]:
    if c in df_u3r.columns:
        df_u3r[c] = df_u3r[c].astype("category")

u3r_meta = {
    "feature_classes": sorted(
        df_u3r["feature_class"].dropna().astype(str).unique().tolist()
    ),
    "features": sorted(df_u3r["feature"].dropna().astype(str).unique().tolist()),
}

df_u3r.to_parquet(OUTDIR / "ltr_u3r_u5.parquet", index=False)
(OUTDIR / "ltr_u3r_u5_meta.json").write_text(json.dumps(u3r_meta, indent=2))






# =========================
# Domains aggregation
# =========================
df = df_dom_raw.copy()

# Basic guards
if "Locus" not in df.columns:
    raise ValueError("Missing 'Locus' in domains TSV")
if "Subfamily" not in df.columns:
    raise ValueError("Missing 'Subfamily' in domains TSV")
if "Type" not in df.columns:
    df["Type"] = pd.NA

# Clean odd whitespaces in headers
df.columns = pd.Index(df.columns).str.replace("\u00A0", " ", regex=False).str.strip()

# Re-key Domains: coord → enumerated (fallback to normalized coord)
df["Locus_norm"] = df["Locus"].map(normalize_internal_name)
df = df.merge(int_map[["Locus_norm", "enum_id"]], on="Locus_norm", how="left")
df["Locus"] = df["enum_id"].where(df["enum_id"].notna(), df["Locus_norm"])

# ---- Coverage detection (robust to naming) ----
norm_map = {re.sub(r'[^a-z]', '', c.lower()): c for c in df.columns}
wanted = [
    "coverageconservedness","coverage","hmmcoverage","domaincoverage",
    "maxcoverage","meancoverage","cov","hmmcov","hmm_coverage","hmmcovpercent",
]
cov_src = None
for key in wanted:
    if key in norm_map:
        cov_src = norm_map[key]
        break

if cov_src is None:
    df["CoverageConservedness"] = pd.NA
else:
    cov = pd.to_numeric(df[cov_src], errors="coerce")
    if cov.dropna().quantile(0.9) > 1.5:  # probably in %
        cov = cov / 100.0
    df["CoverageConservedness"] = cov.clip(lower=0.0, upper=1.0)

# ---- Coarse DomainClass ----
def coarse_domain_class(x: str) -> str:
    if pd.isna(x): return "OTHER"
    u = str(x).upper()
    for k in ["ENV", "POL", "GAG", "ACCESSORY"]:
        if k in u: return k
    m = re.match(r"([A-Z]+)", u)
    return m.group(1) if m else "OTHER"

df["DomainType"] = df["DomainType"].astype(str)
df["DomainClass"] = df["DomainType"].apply(coarse_domain_class)

# ---- Aggregate per Locus × Subfamily × Type ----
agg = (
    df.groupby(["Locus", "Subfamily", "Type"], dropna=False)
      .agg(
          Domains_present=("DomainClass", lambda s: ",".join(sorted(set(s)))),
          N_domain_hits=("DomainClass", "size"),
          MaxCoverage=("CoverageConservedness", "max"),
          MeanCoverage=("CoverageConservedness", "mean"),
      )
      .reset_index()
      .sort_values(["Subfamily", "Locus"])
      .reset_index(drop=True)
)

# ---- Build per-type flags and a compact present-types string ----
df["DomainType_norm"] = (
    df["DomainType"].astype(str).str.upper().str.replace(r"[^A-Z0-9]+", "", regex=True)
      .replace({"GAGCOAT": "GAG"})
)
type_dummies = pd.get_dummies(df["DomainType_norm"])
by_type = type_dummies.groupby(df["Locus"]).max()
by_type.columns = [f"has_{c}" for c in by_type.columns]
domain_types_series = (
    by_type.apply(lambda r: ",".join(sorted([t for t, v in r.items() if v == 1])) or "none", axis=1)
    .rename("DomainTypes_present")
)

# Merge flags into agg
agg = (
    agg.merge(by_type, how="left", left_on="Locus", right_index=True)
       .merge(domain_types_series, how="left", left_on="Locus", right_index=True)
)
for c in by_type.columns:
    agg[c] = agg[c].fillna(0).astype(int)
agg["DomainTypes_present"] = agg["DomainTypes_present"].fillna("none")

# =========================
# Attach INTERNAL metadata (enumerated key preferred)
# =========================
# Decide which columns exist in INTERNAL to use as coords & key
key_col   = "internal_name" if "internal_name" in df_int_full.columns else "name"
start_col = "original_start" if "original_start" in df_int_full.columns else "start"
end_col   = "original_end"   if "original_end"   in df_int_full.columns else "end"

needed_cols = [
    key_col, "RepeatName", "internal_ClassFamily",
    "chrom", start_col, end_col, "strand", "ltr_status"
]
missing = [c for c in needed_cols if c not in df_int_full.columns]
if missing:
    raise ValueError(f"INTERNAL_fully_annotated.tsv is missing expected columns: {missing}")

df_int_min = df_int_full[needed_cols].copy()
df_int_min.rename(columns={
    key_col: "Locus",  # key: enumerated or coordinate-style
    "RepeatName": "Subfamily",
    "internal_ClassFamily": "InternalClass",
    "chrom": "Chrom",
    start_col: "Start",
    end_col: "End",
    "strand": "Strand"
}, inplace=True)


merged = pd.merge(
    df_int_min,
    agg[["Locus","Subfamily","Type","Domains_present","N_domain_hits",
         "MaxCoverage","MeanCoverage","DomainTypes_present"] + by_type.columns.tolist()],
    on="Locus", how="left", suffixes=("_int", "_dom")
)

# Prefer domain Subfamily when present, else internal one
merged["Subfamily"] = merged["Subfamily_dom"].combine_first(merged["Subfamily_int"])
merged.drop(columns=["Subfamily_dom", "Subfamily_int"], inplace=True)

# Defaults & hygiene
merged["Domains_present"] = merged["Domains_present"].fillna("none")
merged["N_domain_hits"]   = merged["N_domain_hits"].fillna(0).astype(int)
for c in ("MaxCoverage", "MeanCoverage"):
    if c in merged.columns:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

# Make sure Chrom/Start/End exist; if any NA, try inferring from Locus when possible
coords = merged["Locus"].apply(lambda x: pd.Series(locus_to_interval(x), index=["Chrom", "Start", "End"]))
for col in ["Chrom", "Start", "End"]:
    if col in merged.columns:
        merged[col] = merged[col].where(merged[col].notna(), coords[col])
    else:
        merged[col] = coords[col]

# Clean subfamily labels for display
merged["Subfamily"] = merged["Subfamily"].map(strip_int_suffix)
merged["Subfamily_clean"] = merged["Subfamily"]

# =========================
# Write outputs
# =========================
merged.to_parquet(OUTDIR / "agg.parquet", index=False)

subfamily_options   = sorted(merged["Subfamily_clean"].dropna().unique().tolist()) if "Subfamily_clean" in merged else []
ltr_status_options  = sorted(merged["ltr_status"].dropna().astype(str).unique().tolist()) if "ltr_status" in merged else []
type_present        = ("Type" in merged) and merged["Type"].notna().any()
type_options        = ["all"] + (sorted(merged["Type"].dropna().unique().tolist()) if type_present else [])

opts = {
    "subfamily_options": subfamily_options,
    "ltr_status_options": ltr_status_options,
    "type_options": type_options
}
(OUTDIR / "domains_meta.json").write_text(json.dumps(opts, indent=2))








print("[OK] Wrote:",
      OUTDIR / "ltr.parquet",
      OUTDIR / "ltr_meta.json",
      OUTDIR / "agg.parquet",
      OUTDIR / "domains_meta.json",
      OUTDIR / "ltr_u3r_u5.parquet",
      OUTDIR / "ltr_u3r_u5_meta.json")

