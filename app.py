#!/usr/bin/env python3
import re
from pathlib import Path

import dash
from dash import dcc, html, dash_table, Input, Output, State, ctx
from dash import callback, Input, Output, no_update
import dash_bootstrap_components as dbc
import dash_bio as dashbio
import pandas as pd
import numpy as np
import json
from dash.dash_table.Format import Format, Scheme

import duckdb
from flask_caching import Cache


## HELPER
def fmt_coord(chrom, start, end, strand=None):
    if pd.isna(chrom) or pd.isna(start) or pd.isna(end):
        return ""
    try:
        s = f"{str(chrom)}:{int(start)}-{int(end)}"
    except Exception:
        return ""
    if strand in ("+","-"):
        s += f" ({strand})"
    return s



PRE = Path("assets/precomputed")

# LTR table + meta
df_ltr = pd.read_parquet(PRE / "ltr.parquet")
ltr_meta = json.loads((PRE / "ltr_meta.json").read_text())

min_len  = int(ltr_meta["min_len"])
max_len  = int(ltr_meta["max_len"])
min_dist = int(ltr_meta["min_dist"])
max_dist = int(ltr_meta["max_dist"])

# Domains (aggregated loci) + meta
agg = pd.read_parquet(PRE / "agg.parquet")
dom_meta = json.loads((PRE / "domains_meta.json").read_text())

# U3/R/U5 + PBS/PPT + signals table + meta
U3R_PARQUET = str(PRE / "ltr_u3r_u5.parquet")
u3r_meta = json.loads((PRE / "ltr_u3r_u5_meta.json").read_text())

u3r_feature_class_options = [
    {"label": s, "value": s} for s in u3r_meta.get("feature_classes", [])
]
u3r_feature_options = [
    {"label": s, "value": s} for s in u3r_meta.get("features", [])
]

# One DuckDB connection reused for on-demand queries
u3r_con = duckdb.connect(database=":memory:")


# Sidebar dropdown options from meta
subfamily_options = [{"label": s, "value": s} for s in dom_meta["subfamily_options"]]
ltr_status_options = [{"label": s, "value": s} for s in dom_meta["ltr_status_options"]]
type_options = [{"label": "All", "value": "all"}] + [
    {"label": t, "value": t} for t in dom_meta["type_options"] if t != "all"
]
# Build "Contains Domain" options from agg's has_* flags
domain_flags = sorted(
    {c.replace("has_", "") for c in agg.columns if c.startswith("has_")}
)
domain_options = [{"label": d, "value": d} for d in domain_flags]

# Add the ACCESSORY meta-bucket if any accessory subtypes exist
if any(d in domain_flags for d in ("ORFX", "ORFQ", "VAP")):
    domain_options = [{"label": "ACCESSORY", "value": "ACCESSORY"}] + domain_options


# --- shared DataTable styles (match domains table look) ---
TABLE_STYLE = {"overflowX": "auto", "maxHeight": "600px", "height": "auto", "minHeight": "500px"}

# For LTR table: keep a fixed vertical size so it doesn't collapse to 1 row
LTR_TABLE_STYLE = {
    "overflowX": "auto",
    "overflowY": "auto",
    "height": "500px",   # or "400px" if you prefer
}

HEADER_STYLE = {
    "backgroundColor": "#f8f9fa",
    "fontWeight": "600",
    "border": "1px solid #e9ecef",
}

# text wraps like your domains table
CELL_STYLE = {
    "minWidth": 80,
    "maxWidth": 320,
    "whiteSpace": "normal",
    "fontSize": "14px",
}

# common zebra striping (apply where you have style_data_conditional)
ZEBRA_STRIPES = [{"if": {"row_index": "odd"}, "backgroundColor": "rgba(0,0,0,0.02)"}]

# clamp Locus column width:
#   - domains table uses column id "Locus"
#   - LTR table uses column id "sequence_name"
LOCUS_CELL_COND = [
    {"if": {"column_id": "Locus"},         "maxWidth": 180, "textAlign": "left"},
    {"if": {"column_id": "Locus_display"}, "maxWidth": 180, "textAlign": "left"},
    {"if": {"column_id": "sequence_name"}, "maxWidth": 180, "textAlign": "left"},
]


# Base UCSC path for GTEx selected tissues bigwig (hg38)
_UCSC_GTEX_BASE = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gtex/cov"

GTEX_RNASEQ = {
    # CNS
    "brain_cortex":      f"{_UCSC_GTEX_BASE}/GTEX-UTHO-3026-SM-3GAFB.Brain_Cortex.RNAseq.bw",
    "brain_frontal":     f"{_UCSC_GTEX_BASE}/GTEX-T5JC-0011-R10A-SM-32PM2.Brain_Frontal_Cortex_BA9.RNAseq.bw",
    "brain_cerebellum":  f"{_UCSC_GTEX_BASE}/GTEX-145MH-2926-SM-5Q5D2.Brain_Cerebellum.RNAseq.bw",
    "brain_hippocampus": f"{_UCSC_GTEX_BASE}/GTEX-1HSKV-0011-R1b-SM-CMKH7.Brain_Hippocampus.RNAseq.bw",

    # Immune
    "whole_blood":       f"{_UCSC_GTEX_BASE}/GTEX-1LG7Z-0005-SM-DKPQ6.Whole_Blood.RNAseq.bw",
    "spleen":            f"{_UCSC_GTEX_BASE}/GTEX-14PKU-0526-SM-6871A.Spleen.RNAseq.bw",
    "ebv_lymphocytes":   f"{_UCSC_GTEX_BASE}/GTEX-1122O-0003-SM-5Q5DL.Cells_EBV-transformed_lymphocytes.RNAseq.bw",

    # Reproductive
    "testis":            f"{_UCSC_GTEX_BASE}/GTEX-1JKYN-1026-SM-CGQG4.Testis.RNAseq.bw",
    "ovary":             f"{_UCSC_GTEX_BASE}/GTEX-ZVT2-0326-SM-5E44G.Ovary.RNAseq.bw",
    "uterus":            f"{_UCSC_GTEX_BASE}/GTEX-1MA7W-1526-SM-DHXKS.Uterus.RNAseq.bw",

    # Barrier / epithelial
    "skin_sun":          f"{_UCSC_GTEX_BASE}/GTEX-1C475-1826-SM-73KWA.Skin_Sun_Exposed_Lower_leg.RNAseq.bw",
    "skin_not_sun":      f"{_UCSC_GTEX_BASE}/GTEX-1JN76-0626-SM-CKZOQ.Skin_Not_Sun_Exposed_Suprapubic.RNAseq.bw",
    "lung":              f"{_UCSC_GTEX_BASE}/GTEX-Y5V5-0826-SM-4VBQD.Lung.RNAseq.bw",

    # Metabolic
    "liver":             f"{_UCSC_GTEX_BASE}/GTEX-Y5LM-0426-SM-4VBRO.Liver.RNAseq.bw",
    "pancreas":          f"{_UCSC_GTEX_BASE}/GTEX-1H1E6-0826-SM-9WG83.Pancreas.RNAseq.bw",

    # Mesenchymal
    "fibroblasts":       f"{_UCSC_GTEX_BASE}/GTEX-117XS-0008-SM-5Q5DQ.Cells_Cultured_fibroblasts.RNAseq.bw",
}

GTEX_LABEL = {
    # CNS
    "brain_cortex":      "Brain – Cortex",
    "brain_frontal":     "Brain – Frontal cortex (BA9)",
    "brain_cerebellum":  "Brain – Cerebellum",
    "brain_hippocampus": "Brain – Hippocampus",

    # Immune
    "whole_blood":       "Whole blood",
    "spleen":            "Spleen",
    "ebv_lymphocytes":   "EBV lymphoblastoid cells",

    # Reproductive
    "testis":            "Testis",
    "ovary":             "Ovary",
    "uterus":            "Uterus",

    # Barrier / epithelial
    "skin_sun":          "Skin – Sun-exposed (leg)",
    "skin_not_sun":      "Skin – Not sun-exposed",
    "lung":              "Lung",

    # Metabolic
    "liver":             "Liver",
    "pancreas":          "Pancreas",

    # Mesenchymal
    "fibroblasts":       "Cultured fibroblasts",
}




# --------------------
# App (create first so cache can bind to app.server)
# --------------------
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    update_title=None,          # ← stop switching to "Updating..."
)
app.title = "HERVarium"

# Small cache so repeated queries are instant (memory store is fine)
cache = Cache(config={"CACHE_TYPE": "SimpleCache", "CACHE_DEFAULT_TIMEOUT": 60})
cache.init_app(app.server)

# --------------------
# Backend connections & helpers
# --------------------

# Safe global bounds for LTR length

min_len = int(min_len)
max_len = int(max_len)
if min_len >= max_len:
    max_len = min_len + 1

print(f"[HERVarium] LTR length bounds: {min_len}..{max_len}")

# Parse IGV locus strings like "chr10:6,824,465-6,833,272"
LOCI_RE = re.compile(r"^(chr[0-9XYM]+):([\d,]+)-([\d,]+)$")
def parse_locus(s: str):
    if not s:
        return None
    m = LOCI_RE.match(s.strip())
    if not m:
        return None
    chrom = m.group(1)
    start = int(m.group(2).replace(",", ""))
    end   = int(m.group(3).replace(",", ""))
    if end < start:
        start, end = end, start
    return chrom, start, end





# --------------------
# UI components
# --------------------
# Navbar with logo
navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                dbc.Row(
                    [
                        dbc.Col(
                            html.Img(
                                src="/assets/logos/hervarium_logo.png",  # use .png if you add it later
                                className="brand-logo",
                                alt="HERVarium logo",
                                draggable="false",
                            ),
                            width="auto",
                        ),
                        dbc.Col(
                            dbc.NavbarBrand(
                                "HERVarium — Genome-wide atlas of HERV domains, LTR motifs, and U3–R–U5 architecture",
                                className="ms-1 brand-text",
                            ),
                            width="auto",
                        ),
                    ],
                    align="center",
                    className="brand-wrap",
                ),
                href="#",
                style={"textDecoration": "none"},
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
        ],
        fluid=True,
    ),
    color="dark",
    dark=True,
    className="mb-3 shadow-sm",
)

# --------------------
# Sidebar — tidy, card-based
# --------------------
sidebar = dbc.Col(
    dbc.Card(
        dbc.CardBody(
            [
                html.H5("Chromatin accessibility filters", className="card-title"),

                # DNase cell-type (ENCODE/UW) selector
                # DNase cell-type (ENCODE/UW) selector + GTEx RNA-seq tissue
                html.Div(
                    [
                        dbc.Label("DNase I HS peaks (ENCODE/UW) — cell type"),
                        dcc.Dropdown(
                            id="sel-dnase",
                            options=[
                                {"label": "GM12878 (B-lymphoblastoid)", "value": "gm12878"},
                                {"label": "K562 (erythroleukemia)",      "value": "k562"},
                                {"label": "HepG2 (liver)",               "value": "hepg2"},
                                {"label": "H7-hESC",                     "value": "h7hesc"},
                                {"label": "HeLa-S3 (cervical cancer)",   "value": "helas3"},
                                {"label": "CD14+ monocytes",             "value": "monocytescd14ro01746"},
                                {"label": "CD20+ B cells",               "value": "cd20ro01778"},
                                {"label": "Jurkat (T cells)",            "value": "jurkat"},
                                {"label": "Th1",                         "value": "th1"},
                                {"label": "Th2",                         "value": "th2"},
                                {"label": "Be2C (neuroblastoma)",        "value": "be2c"},
                                {"label": "M059J (glioblastoma)",        "value": "m059j"},
                                {"label": "HUVEC (endothelial)",         "value": "huvec"},
                                {"label": "HSMM (muscle differentiation)","value": "hsmm"},
                                {"label": "Wi-38 (fibroblast)",          "value": "wi38"},
                            ],
                            value="gm12878",
                            clearable=False,
                        ),

                        dbc.Checkbox(
                            id="chk-show-dnase",
                            value=False,
                            label="Show DNase hotspots track",
                            className="mt-2",
                        ),

                        html.Hr(className="my-2"),

                        dbc.Label("GTEx RNA-seq — tissue"),
                        dcc.Dropdown(
                            id="sel-gtex",
                            options=[
                                {"label": label, "value": key}
                                for key, label in GTEX_LABEL.items()
                            ],
                            value="whole_blood",   # or "brain_cortex" if you prefer CNS as default
                            clearable=False,
                        ),
                        dbc.Checkbox(
                            id="chk-show-gtex",
                            value=False,
                            label="Show GTEx RNA-seq track",
                            className="mt-1",
                        ),
                    ],
                    className="mb-3",
                ),

                
                html.Hr(className="my-3"),
                html.H5("Domain filters", className="card-title"),

                # Domain: Subfamily
                html.Div(
                    [
                        dbc.Label("Subfamily"),
                        dcc.Dropdown(
                            id="flt-subfamily",
                            options=subfamily_options,
                            placeholder="Select subfamily",
                            multi=True,
                            persistence=True,
                            persistence_type="session",
                        ),
                    ],
                    className="mb-3",
                ),

                # Domain: Contains Domain
                html.Div(
                    [
                        dbc.Label("Contains Domain"),
                        dcc.Dropdown(
                            id="flt-domain",
                            options=domain_options,
                            placeholder="Select domain class(es)",
                            multi=True,
                            persistence=True,
                            persistence_type="session",
                        ),
                    ],
                    className="mb-3",
                ),

                # Domain: Min Coverage
                html.Div(
                    [
                        dbc.Label("Min Coverage (HMM coverage)"),
                        dcc.Slider(
                            id="flt-mincov",
                            min=0.0, max=1.0, step=0.01, value=0.0,
                            marks={0: "0", 0.5: "0.5", 1: "1"},
                            tooltip={"placement": "bottom"},
                        ),
                    ],
                    className="mb-2",
                ),

                # At least one domain (N_domain_hits ≥ 1)
                html.Div(
                    [
                        dbc.Switch(
                            id="flt-has-domain",
                            value=False,
                            label="At least one domain (N ≥ 1)",
                            persistence=True, persistence_type="session",
                        ),
                    ],
                    className="mb-3",
                ),

                # LTR status (multi)
                html.Div(
                    [
                        dbc.Label("LTR status"),
                        dcc.Dropdown(
                            id="flt-ltrstatus",
                            options=ltr_status_options,   # built above from agg["ltr_status"]
                            multi=True,
                            placeholder="Any",
                            persistence=True, persistence_type="session",
                        ),
                    ],
                    className="mb-3",
                ),

                html.Hr(className="my-3"),
                html.H5("LTR filters", className="card-title"),


                # LTR Subfamily (multi)
                html.Div([
                    dbc.Label("LTR Subfamily"),
                    dcc.Dropdown(
                        id="flt-ltr-subfamily",
                        options=[{"label": s, "value": s} for s in sorted(df_ltr["subfamily"].dropna().unique().tolist())],
                        multi=True, placeholder="Select LTR subfamily",
                        persistence=True, persistence_type="session",
                    ),
                ], className="mb-3"),

                # Type: solo_LTR / 5prime_LTR / 3prime_LTR
                html.Div([
                    dbc.Label("Type of LTR"),
                    dcc.Dropdown(
                        id="flt-ltr-type",
                        options=[{"label": s, "value": s} for s in sorted(df_ltr["type"].dropna().unique().tolist())],
                        multi=True, placeholder="Select type",
                        persistence=True, persistence_type="session",
                    ),
                ], className="mb-3"),

                # Internal region subfamily
                html.Div([
                    dbc.Label("Internal region subfamily"),
                    dcc.Dropdown(
                        id="flt-ltr-intsub",
                        options=[{"label": s, "value": s} for s in sorted(df_ltr["int_subfamily"].dropna().unique().tolist())],
                        multi=True, placeholder="Select internal subfamily",
                        persistence=True, persistence_type="session",
                    ),
                ], className="mb-3"),

                # Distance to TSS
                html.Div([
                    dbc.Label("Distance to nearest TSS (bp)"),
                    dbc.InputGroup([
                        dbc.InputGroupText("min"),
                        dbc.Input(id="flt-ltr-dist-min", type="number", step=1, value=min_dist, debounce=True),
                        dbc.InputGroupText("max"),
                        dbc.Input(id="flt-ltr-dist-max", type="number", step=1, value=max_dist, debounce=True),
                    ]),
                    html.Small("Absolute genomic distance reported in the table.", className="text-muted"),
                ], className="mb-3"),

                # LTR length
                html.Div([
                    dbc.Label("LTR length (bp)"),
                    dbc.InputGroup([
                        dbc.InputGroupText("min"),
                        dbc.Input(id="flt-ltr-len-min", type="number", step=1, value=min_len, debounce=True),
                        dbc.InputGroupText("max"),
                        dbc.Input(id="flt-ltr-len-max", type="number", step=1, value=max_len, debounce=True),
                    ]),
                    html.Small("Filter by genomic length: end − start + 1", className="text-muted"),
                ], className="mb-3"),

                # Minimum #motifs
                html.Div([
                    dbc.Label("Min #Motifs"),
                    dbc.Input(id="flt-ltr-nmotifs", type="number", min=0, step=1, value=0, debounce=True),
                ], className="mb-3"),

                html.Hr(className="my-3"),
                html.H5("U3/R/U5 feature filters", className="card-title"),

                # U3/R/U5 feature class
                html.Div(
                    [
                        dbc.Label("Feature class"),
                        dcc.Dropdown(
                            id="flt-u3r-class",
                            options=u3r_feature_class_options,
                            multi=True,
                            placeholder="SEGMENT / PBS_PPT / SIGNAL",
                        ),
                    ],
                    className="mb-3",
                ),

                # U3/R/U5 feature
                html.Div(
                    [
                        dbc.Label("Feature"),
                        dcc.Dropdown(
                            id="flt-u3r-feature",
                            options=u3r_feature_options,
                            multi=True,
                            placeholder="U3, R, U5, PBS, PPT, TSS, ...",
                        ),
                    ],
                    className="mb-3",
                ),

                # U3/R/U5 score + LOW_CONF toggle
                html.Div(
                    [
                        dbc.Label("Min score"),
                        dbc.Input(
                            id="flt-u3r-min-score",
                            type="number",
                            step=0.1,
                            min=0,
                            debounce=True,
                        ),
                        dbc.Checkbox(
                            id="flt-u3r-hide-lowconf",
                            value=False,
                            label="Hide LOW_CONF segments",
                            className="mt-1",
                        ),
                    ],
                    className="mb-3",
                ),

                dbc.Button(
                    "Clear filters",
                    id="btn-clear",
                    color="secondary",
                    size="sm",
                    className="mt-1",
                ),


            ]
        ),
        className="shadow-sm",
    ),
    width=3,
)


# IGV tracks
# Base UCSC path for ENCODE UW DNase hotspot bigBeds (hg38)
_UCSC_DNASE_BASE = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/wgEncodeRegDnase"

DNASE_BIGBED = {
    "gm12878": f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwGm12878Hotspot.broadPeak.bb",
    "k562":    f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwK562Hotspot.broadPeak.bb",
    "hepg2":   f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwHepg2Hotspot.broadPeak.bb",
    "h7hesc":  f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwH7hescHotspot.broadPeak.bb",
    "helas3":  f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwHelas3Hotspot.broadPeak.bb",
    "monocytescd14ro01746": f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwMonocytescd14ro01746Hotspot.broadPeak.bb",
    "cd20ro01778":          f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwCd20ro01778Hotspot.broadPeak.bb",
    "jurkat":  f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwJurkatHotspot.broadPeak.bb",
    "th1":     f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwTh1Hotspot.broadPeak.bb",
    "th2":     f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwTh2Hotspot.broadPeak.bb",
    "be2c":    f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwBe2cHotspot.broadPeak.bb",
    "m059j":   f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwM059jHotspot.broadPeak.bb",
    "huvec":   f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwHuvecHotspot.broadPeak.bb",
    "hsmm":    f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwHsmmHotspot.broadPeak.bb",
    "wi38":    f"{_UCSC_DNASE_BASE}/wgEncodeRegDnaseUwWi38Hotspot.broadPeak.bb",
}

# (Optional) pretty names to use in the IGV track title
DNASE_LABEL = {
    "gm12878": "GM12878 (B-lymphoblastoid)",
    "k562":    "K562 (erythroleukemia)",
    "hepg2":   "HepG2 (liver)",
    "h7hesc":  "H7-hESC",
    "helas3":  "HeLa-S3 (cervical cancer)",
    "monocytescd14ro01746": "CD14+ monocytes",
    "cd20ro01778":          "CD20+ B cells",
    "jurkat":  "Jurkat (T cells)",
    "th1":     "Th1",
    "th2":     "Th2",
    "be2c":    "Be2C (neuroblastoma)",
    "m059j":   "M059J (glioblastoma)",
    "huvec":   "HUVEC (endothelial)",
    "hsmm":    "HSMM (muscle differentiation)",
    "wi38":    "Wi-38 (fibroblast)",
}


BASE_TRACKS = [
    {
        "name": "GENCODE genes (symbols)",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/gencode/gencode.v48.genesymbols.bb",
        "displayMode": "EXPANDED",
        "height": 150,
        "color": "#7F7F7F",   # neutral gray
    },
    {
        "name": "HERV internal regions (GyDB)",
        "type": "annotation",
        "format": "bed",
        "url": "/assets/internals/HERV_internal_simplified.bed",
        "displayMode": "EXPANDED",
        "visibilityWindow": 200000,
        "height": 50,
        "color": "#0072B2",   # deep blue (pair with domains)
    },
    {
        "name": "HERV domains (BED)",
        "type": "annotation",
        "format": "bed",
        "url": "/assets/internals/HERV_internal_domains_simplified.bed",
        "displayMode": "EXPANDED",
        "visibilityWindow": 200000,
        "height": 100,
        "useScore": True,
        "color": "#56B4E9",   # sky blue (lighter of the pair)
    },
    {
        "name": "HERV LTRs",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/ltr/ERV_ltr_merged.simplified.bb",
        "displayMode": "EXPANDED",
        "visibilityWindow": 200000,
        "height": 100,
        "color": "#D55E00",   # vermillion (pair with LTR annotations)
    },

    {
        "name": "LTR U3/R/U5 segments",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/ltr/segments/HERV_LTR_U3_R_U5_segments_allconf.bb",  # point to the .bb
        "displayMode": "EXPANDED",
        "visibilityWindow": 200000,   # narrower window → lighter per-view load
        "height": 90,
        "color": "#CC79A7",   # reddish-purple
    },

    {
        "name": "LTR PBS/PPT",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/ltr/segments/HERV_LTR_U3_R_U5_PBS_PPT.bb",
        "displayMode": "EXPANDED",
        "visibilityWindow": 200000,
        "height": 60,
        "color": "#009E73",   # teal
    },

    {
        "name": "LTR promoter/PAS signals",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/ltr/segments/HERV_LTR_U3_R_U5_signals.bb", 
        "displayMode": "EXPANDED",
        "visibilityWindow": 50000,   
        "height": 100,
        "color": "#F0E442",
    },

    {
        "name": "LTR TFBMs",
        "type": "annotation",
        "format": "bigbed",
        "url": "/assets/ltr/tfbm/fimo_parsed.bb",
        "displayMode": "SQUISHED",
        "visibilityWindow": 200000,
        "height": 300,
        "color": "#E69F00",   # orange (lighter warm of the pair)
    },
]



def dnase_track(cellkey: str):
    url = DNASE_BIGBED.get(cellkey)
    if not url:
        return None
    title = DNASE_LABEL.get(cellkey, f"DNase hotspots ({cellkey})")
    return {
        "name": f"ENCODE DNase hotspots — {title}",
        "type": "annotation",
        "format": "bigbed",
        "url": url,
        "displayMode": "SQUISHED",
        "height": 70,
        "visibilityWindow": 500000,
        "color": "#8A5E9E",  # paired with GTEx yellow
        "alpha": 0.95,       # optional: softens the fill a touch
    }

def gtex_track(tissue_key: str):
    """
    Build an IGV bigWig track for GTEx RNA-seq coverage in a selected tissue.
    """
    url = GTEX_RNASEQ.get(tissue_key)
    if not url:
        return None

    title = GTEX_LABEL.get(tissue_key, tissue_key)
    return {
        "name": f"GTEx RNA-seq — {title}",
        "type": "wig",        # IGV bigWig
        "format": "bigwig",
        "url": url,
        "autoscale": True,
        "height": 80,
        "color": "#EBCFE4",   # soft lilac (paired with DNase purple)
    }


# IGV browser
igv_browser = dashbio.Igv(
    id="genome-browser",
    reference={
        "id": "hg38",
        "name": "Human (GRCh38)",
        "fastaURL": "/assets/genome/GRCh38.primary_assembly.genome.fa",
        "indexURL": "/assets/genome/GRCh38.primary_assembly.genome.fa.fai",
    },
    locus="chr10:6,823,465-6,834,272",
    minimumBases=1,
    tracks=BASE_TRACKS,  # no DNase or GTEx at startup
)


# Main table columns
table_columns = [
    {"name": "Locus", "id": "Locus_display"},  # visible: "chr2:330-560"
    {"name": "Locus (ID)", "id": "Locus"},     # hidden (used for linking)
    {"name": "Subfamily", "id": "Subfamily"},
    {"name": "LTR status", "id": "ltr_status"},
    {"name": "Domains_present", "id": "Domains_present"},
    {"name": "N_domain_hits", "id": "N_domain_hits", "type": "numeric"},
    {"name": "MaxCoverage", "id": "MaxCoverage", "type": "numeric",
     "format": Format(precision=3, scheme=Scheme.fixed)},
    {"name": "MeanCoverage", "id": "MeanCoverage", "type": "numeric",
     "format": Format(precision=3, scheme=Scheme.fixed)},
    {"name": "Chrom", "id": "Chrom"},
    {"name": "Start", "id": "Start", "type": "numeric"},
    {"name": "End", "id": "End", "type": "numeric"},
]



motif_table = [
    {"name": "Locus",           "id": "Locus_display"},  # visible coords
    {"name": "sequence_name",   "id": "sequence_name"},  # optional: keep but hide (below)
    {"name": "Subfamily",       "id": "subfamily"},
    {"name": "ClassFamily",     "id": "ltr_ClassFamily"},
    {"name": "Chrom",           "id": "chrom"},
    {"name": "Start",           "id": "ltr_start", "type": "numeric"},
    {"name": "End",             "id": "ltr_end",   "type": "numeric"},
    {"name": "Length",          "id": "ltr_len",   "type": "numeric"},
    {"name": "Strand",          "id": "ltr_strand"},
    {"name": "#Motifs",         "id": "n_motifs",  "type": "numeric"},
    {"name": "Type",            "id": "type"},
    {"name": "Internal name",   "id": "internal_name"},   # enumerated ID used for linking
    {"name": "Internal subfam", "id": "int_subfamily"},
    {"name": "dist_to_tss",     "id": "dist_to_tss", "type": "numeric"},
    {"name": "nearest_gene",    "id": "nearest_gene"},
    {"name": "tss_bin",         "id": "tss_bin"},
]


u3r_table_columns = [
    {"name": "Feature",         "id": "feature"},
    {"name": "Class",           "id": "feature_class"},
    {"name": "Feature locus",   "id": "feature_locus"},
    {"name": "Score",           "id": "score", "type": "numeric"},
    {"name": "Segment conf.",   "id": "segment_conf"},
    {"name": "LTR sequence",    "id": "sequence_name"},
    {"name": "LTR subfamily",   "id": "subfamily"},
    {"name": "LTR type",        "id": "type"},
    {"name": "Internal subfam", "id": "int_subfamily"},
    {"name": "Nearest gene",    "id": "nearest_gene"},
    {"name": "dist_to_tss",     "id": "dist_to_tss", "type": "numeric"},
]



# ------------------------
# Content column (main area)
# ------------------------
content = dbc.Col(
    [
        # Locus Browser
        dbc.Card(
            dbc.CardBody(
                [
                    html.H4("Locus Browser", className="mb-3"),
                    igv_browser,
                ]
            ),
            className="mb-3 shadow-sm",
        ),

        # Annotated HERV Loci (Domains)
        dbc.Card(
            dbc.CardBody(
                [

                    html.H5("Annotated HERV Loci", className="mb-2"),
                    dash_table.DataTable(
                        id="tbl-loci",
                        data=agg.to_dict("records"),
                        columns=table_columns,
                        hidden_columns=["Locus"],
                        page_size=15,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        row_selectable="single",
                        selected_rows=[],
                        export_format="csv",
                        export_headers="display",
                        style_table=TABLE_STYLE,
                        style_cell=CELL_STYLE,
                        style_header=HEADER_STYLE,
                        style_data_conditional=ZEBRA_STRIPES + [
                            {"if": {"column_id": "MaxCoverage"}, "textAlign": "right"},
                            {"if": {"column_id": "MeanCoverage"}, "textAlign": "right"},
                            {"if": {"column_id": "N_domain_hits"}, "textAlign": "right"},
                        ],
                        fixed_rows={"headers": True},
                        virtualization=True,
                        style_cell_conditional=LOCUS_CELL_COND,
                        tooltip_data=[
                            {"Locus": {"value": str(row["Locus"]), "type": "text"}}
                            for row in agg.to_dict("records")
                        ],
                        tooltip_duration=None,
                    ),
                    # --- link domains → LTRs controls ---
                    html.Div(
                        [
                            dbc.Button(
                                "Link to LTRs (apply current domain filters/selection)",
                                id="btn-link-ltrs",
                                color="primary",
                                size="sm",
                                className="me-2",
                            ),
                            dbc.Checkbox(
                                id="chk-link-use-selected",
                                value=True,
                                label="Use selected rows (otherwise, use all filtered rows)",
                                className="me-3",
                            ),
                            html.Span(id="txt-link-summary", className="text-muted"),
                        ],
                        className="mb-2",
                    ),
                    dbc.Button(
                        "Unlink",
                        id="btn-unlink",
                        color="secondary",
                        size="sm",
                        outline=True,
                        className="ms-2",
                    ),

                ]
            ),
            className="mb-3 shadow-sm",
        ),

        # LTR Elements Table
        dbc.Card(
            dbc.CardBody(
                [
                    html.H5("Annotated LTR Elements", className="mb-2"),
                    dash_table.DataTable(
                        id="ltr-table",
                        columns=motif_table,
                        hidden_columns=["sequence_name"],
                        data=[],
                        page_size=15,
                        page_current=0,
                        page_action="custom",
                        page_count=1,
                        sort_action="custom",
                        sort_mode="multi",
                        sort_by=[],
                        export_format="csv",
                        export_headers="display",
                        filter_action="none",
                        row_selectable="single",
                        selected_rows=[],
                        style_table=LTR_TABLE_STYLE,
                        style_cell=CELL_STYLE,
                        style_header=HEADER_STYLE,
                        style_data_conditional=ZEBRA_STRIPES + [
                            {"if": {"column_id": "ltr_start"}, "textAlign": "right"},
                            {"if": {"column_id": "ltr_end"},   "textAlign": "right"},
                            {"if": {"column_id": "ltr_len"},   "textAlign": "right"},
                            {"if": {"column_id": "n_motifs"},  "textAlign": "right"},
                        ],
                        fixed_rows={"headers": True},
                        virtualization=False,
                        style_cell_conditional=LOCUS_CELL_COND,
                        tooltip_data=[],
                        tooltip_duration=None,
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Link to U3/R/U5, PBS/PPT (apply current LTR filters/selection)",
                                id="btn-link-u3r",
                                color="primary",
                                size="sm",
                                className="me-2 mt-2",
                            ),
                            dbc.Checkbox(
                                id="chk-link-u3r-use-selected",
                                value=True,
                                label="Use selected LTRs (otherwise, all filtered)",
                                className="me-3 mt-2",
                            ),
                            dbc.Button(
                                "Unlink U3/R/U5 table",
                                id="btn-unlink-u3r",
                                color="secondary",
                                size="sm",
                                outline=True,
                                className="mt-2 me-2",
                            ),
                            html.Span(
                                id="txt-link-u3r-summary",
                                className="text-muted ms-1 mt-2",
                            ),
                        ],
                        className="mb-2",
                    ),
                ]
            ),
            className="mb-3 shadow-sm",
        ),

        # U3/R/U5 + PBS/PPT + promoter/PAS annotations
        dbc.Card(
            dbc.CardBody(
                [
                    html.H5(
                        "U3/R/U5 segments, PBS/PPT and promoter/PAS signals",
                        className="mb-2",
                    ),
                    dash_table.DataTable(
                        id="u3r-table",
                        columns=u3r_table_columns,
                        data=[],
                        page_size=15,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        export_format="csv",
                        export_headers="display",
                        style_table=TABLE_STYLE,
                        style_cell=CELL_STYLE,
                        style_header=HEADER_STYLE,
                        style_data_conditional=ZEBRA_STRIPES,
                        fixed_rows={"headers": True},
                        virtualization=False,
                    ),
                ]
            ),
            className="mb-3 shadow-sm",
        ),



        # Footer
        html.Div(
            [
                html.Span("Data availability:\u00A0", className="fw-bold"),
                html.A("LTR regulatory atlas\u00A0", 
                    href="https://doi.org/10.5281/zenodo.17602210",
                    target="_blank"),
                html.Span(" · "),
                html.A("\u00A0Internal domain annotation", 
                    href="https://doi.org/10.5281/zenodo.16318927",
                    target="_blank"),
            ],
            style={"display": "flex", "justifyContent": "center", "marginTop": "5px"},
            className="text-muted small mt-3 mb-4",
        ),
        html.Hr(style={"width": "60%", "margin": "10px auto", "borderTop": "1px solid #ccc"}),
        html.Div(
            [
                html.Div("© HERVarium — Genome-wide atlas of HERV domains, LTR motifs, and U3–R–U5 architecture", style={"marginTop": "15px", "display": "flex", "justifyContent": "center"},), 
                html.Div("Developed by the Functional Genomics Team, CNAG", style={"display": "flex", "justifyContent": "center"},),
                html.Div("AGAUR-FI predoctoral grant (2025 FI-1 00642, Joan Oró) and the European Social Fund Plus.", style={"display": "flex", "justifyContent": "center"},),
                html.Div("CNAG, Spanish Ministry of Science and Innovation (PI19/01772), ERDF, Generalitat de Catalunya.", style={"display": "flex", "justifyContent": "center"},),
                # Logos row
                html.Div(
                    [
                        html.Img(src="/assets/logos/logo_cnag.jpg", height="50px", style={"marginRight": "15px"}),
                        html.Img(src="/assets/logos/logo_generalitat.png", height="50px"),
                        html.Img(src="/assets/logos/logo_eu.png", height="50px", style={"marginRight": "15px"}),
                    ],
                    style={"marginTop": "15px", "display": "flex", "justifyContent": "center"},
                ),
            ],
            className="text-muted small mt-3 mb-4",
        ), 

    ],
    
    width=9,
)



# Full layout
app.layout = html.Div(
    [
        navbar,
        dbc.Container(dbc.Row([sidebar, content]), fluid=True),

        # --- helpers for IGV viewport preservation ---
        dcc.Store(id="store-locus", data=None),
        dcc.Store(id="store-linked-internal-names", data=[]),
        dcc.Store(id="store-linked-ltr-seqnames", data=[]), 
    ]
)


# --------------------
# Helpers & Callbacks
# --------------------
from dash.exceptions import PreventUpdate
import urllib.parse
import pandas as pd
import re

import re
import pandas as pd
import numpy as np

_AND_SPLIT_RE = re.compile(r"\s*&&\s*")
_PARENS_RE = re.compile(r"^\((.*)\)$")  # strip one outer (...) pair


def _strip_outer_parens(s: str) -> str:
    s = (s or "").strip()
    while True:
        m = _PARENS_RE.match(s)
        if not m:
            return s
        s = m.group(1).strip()


def _apply_bar_filter(df: pd.DataFrame, filter_query: str) -> pd.DataFrame:
    """
    Apply Dash DataTable 'filter_query' (custom filtering) to a pandas DataFrame.
    Supports AND (&&). Ignores OR (||) for now.
    """
    if not filter_query:
        return df

    # Dash may include parentheses and variable spacing
    parts = [_strip_outer_parens(p) for p in _AND_SPLIT_RE.split(filter_query) if p.strip()]

    for part in parts:
        col, op, val = _split_filter_part(part)  # <-- use your existing helper
        if not col or op is None:
            continue
        if col not in df.columns:
            continue

        s = df[col]

        # contains / icontains
        if op in ("contains", "icontains"):
            text = s.astype("string")
            mask = text.str.contains(val, case=(op != "icontains"), na=False)

        # numeric comparisons (try to coerce even if dtype is object/Int64)
        elif op in (">=", ">", "<=", "<", "=", "!="):
            # try numeric first
            num = pd.to_numeric(s, errors="coerce")
            v = pd.to_numeric(pd.Series([val]), errors="coerce").iloc[0]

            if not pd.isna(v) and num.notna().any():
                if op == "=":
                    mask = num == v
                elif op == "!=":
                    mask = num != v
                elif op == ">":
                    mask = num > v
                elif op == ">=":
                    mask = num >= v
                elif op == "<":
                    mask = num < v
                elif op == "<=":
                    mask = num <= v
            else:
                # fallback to string equality/inequality only
                text = s.astype("string")
                if op == "=":
                    mask = text == val
                elif op == "!=":
                    mask = text != val
                else:
                    # can't do >,< on non-numeric strings
                    continue
        else:
            continue

        df = df[mask]

    return df




def apply_filters(df_in: pd.DataFrame, subfamilies, domains, mincov,
                  must_have_domain=False, ltr_statuses=None):
    """Apply filters for the HERV domains table (agg DataFrame)."""
    out = df_in.copy()

    # Subfamily
    if subfamilies:
        out = out[out["Subfamily_clean"].isin(subfamilies)]

    # Domain TYPE filter
    if domains:
        for t in (domains or []):
            tt = str(t).upper()
            if tt == "ACCESSORY":
                # any of ORFX/ORFQ/VAP
                acc_cols = [c for c in ["has_ORFX","has_ORFQ","has_VAP"] if c in out.columns]
                if acc_cols:
                    out = out[out[acc_cols].sum(axis=1) >= 1]
                else:
                    # fallback to string if flags missing (shouldn't happen)
                    out = out[out["DomainTypes_present"].str.contains("ORFX|ORFQ|VAP", na=False)]
            else:
                col = f"has_{tt}"
                if col in out.columns:
                    out = out[out[col] == 1]
                else:
                    # fallback to string matching of DomainTypes_present
                    out = out[out["DomainTypes_present"].str.contains(
                        fr"(?:^|,){re.escape(tt)}(?:,|$)", na=False)]


    # Min coverage (treat NaN as 0)
    if "MaxCoverage" in out.columns and mincov is not None:
        out = out[(out["MaxCoverage"].fillna(0) >= float(mincov))]

    # NEW: must have at least one domain
    if must_have_domain:
        out = out[out["N_domain_hits"] >= 1]

    # NEW: LTR status filter
    if ltr_statuses:
        out = out[out["ltr_status"].isin(ltr_statuses)]

    return out


# ---------- Domains table ----------
@app.callback(
    Output("tbl-loci", "data"),
    Input("flt-subfamily", "value"),
    Input("flt-domain", "value"),
    Input("flt-mincov", "value"),
    Input("flt-has-domain", "value"),   # NEW
    Input("flt-ltrstatus", "value"),    # NEW
)
def update_table(subfamilies, domains, mincov, must_have_domain, ltr_statuses):
    subfamilies = subfamilies or []
    domains = domains or []
    ltr_statuses = ltr_statuses or []
    filtered = apply_filters(agg, subfamilies, domains, mincov,
                             must_have_domain=bool(must_have_domain),
                             ltr_statuses=ltr_statuses)

    # Build a pretty coord string for display; fall back to ID if coords are missing
    def _fmt(row):
        c, s, e = row.get("Chrom"), row.get("Start"), row.get("End")
        if pd.isna(c) or pd.isna(s) or pd.isna(e):
            return str(row.get("Locus", ""))
        try:
            return f"{str(c)}:{int(s)}-{int(e)}"
        except Exception:
            return str(row.get("Locus", ""))

    filtered = filtered.copy()
    # AFTER (with strand)
    filtered["Locus_display"] = filtered.apply(
        lambda r: fmt_coord(r.get("Chrom"), r.get("Start"), r.get("End"), r.get("Strand")),
        axis=1
    )


    return filtered.to_dict("records")




# From Domains table -> jump IGV
@app.callback(
    Output("genome-browser", "locus", allow_duplicate=True),
    Output("store-locus", "data", allow_duplicate=True),   # <-- add allow_duplicate 
    Input("tbl-loci", "selected_rows"),
    State("tbl-loci", "data"),
    prevent_initial_call=True,
)
def jump_igv_loci(selected_rows, table_data):
    from dash.exceptions import PreventUpdate
    import pandas as pd
    if not table_data or not selected_rows:
        raise PreventUpdate
    row = table_data[selected_rows[0]]
    chrom, start, end = row.get("Chrom"), row.get("Start"), row.get("End")
    if not chrom or pd.isna(start) or pd.isna(end):
        raise PreventUpdate
    pad = max(200, int(0.05 * (end - start)))
    locus = f"{chrom}:{max(1, int(start - pad))}-{int(end + pad)}"
    return locus, locus

# Swap DNase track WITHOUT losing the current IGV viewport
from dash import no_update

@app.callback(
    Output("genome-browser", "tracks"),
    Input("sel-dnase", "value"),
    Input("chk-show-dnase", "value"),
    Input("sel-gtex", "value"),
    Input("chk-show-gtex", "value"),
)
def update_igv_tracks(cellkey, show_dnase, gtex_key, show_gtex):
    """
    Build IGV tracks dynamically:
      - Start from BASE_TRACKS (static).
      - Optionally add DNase hotspots for the selected cell type.
      - Optionally add GTEx RNA-seq bigWig for the selected tissue.
    """
    tracks = list(BASE_TRACKS)

    # GTEx RNA-seq coverage
    if show_gtex:
        key = gtex_key or "whole_blood"
        t_gtex = gtex_track(key)
        if t_gtex is not None:
            tracks.append(t_gtex)

    # ENCODE DNase hotspots
    if show_dnase:
        key_dnase = (cellkey or "gm12878").lower()
        t_dnase = dnase_track(key_dnase)
        if t_dnase is not None:
            tracks.append(t_dnase)

    return tracks






# ---------- Server-side LTR table ----------
def _split_filter_part(filter_part: str):
    """Parse a single filter expression from DataTable."""
    if not filter_part or not filter_part.strip():
        return None, None, None
    s = filter_part.strip()
    col = None
    if s.startswith("{") and "}" in s:
        col = s[1:s.index("}")]
        s = s[s.index("}") + 1:].lstrip()
    ops = ["contains", "icontains", ">=", "<=", ">", "<", "!=", "="]
    for op in ops:
        if s.lower().startswith(op):
            val = s[len(op):].strip()
            val = urllib.parse.unquote(val).strip().strip('"').strip("'")
            return col, op, val
    return None, None, None


@callback(
    Output("ltr-table", "data"),
    Output("ltr-table", "page_count"),
    Input("ltr-table", "page_current"),
    Input("ltr-table", "page_size"),
    Input("ltr-table", "sort_by"),
    Input("flt-ltr-subfamily", "value"),
    Input("flt-ltr-len-min", "value"),
    Input("flt-ltr-len-max", "value"),
    Input("flt-ltr-nmotifs", "value"),
    # NEW inputs you added earlier for type/int_subfam/dist:
    Input("flt-ltr-type", "value"),
    Input("flt-ltr-intsub", "value"),
    Input("flt-ltr-dist-min", "value"),
    Input("flt-ltr-dist-max", "value"),
    # <-- NEW: the link payload
    Input("store-linked-internal-names", "data"),
)
def ltr_table_page(page_current, page_size, sort_by,
                   subf_sel, len_min, len_max, min_nmotifs,
                   type_sel, intsub_sel, dist_min, dist_max,
                   linked_internal_names):

    page_current = page_current or 0
    page_size    = page_size or 15

    df = df_ltr  # your in-memory LTR dataframe

    # --- apply link (optional) ---
    if linked_internal_names:
        df = df[df["internal_name"].isin(linked_internal_names)]


    # --- existing filters (subfamily/type/int_subfam/n_motifs/length/dist/etc.) ---
    if subf_sel:
        df = df[df["subfamily"].isin(subf_sel)]
    if type_sel:
        df = df[df["type"].isin(type_sel)]
    if intsub_sel:
        df = df[df["int_subfamily"].isin(intsub_sel)]
    if min_nmotifs is not None and int(min_nmotifs) > 0:
        df = df[df["n_motifs"] >= int(min_nmotifs)]

    # ranges...
    lmin = int(len_min) if len_min is not None else min_len
    lmax = int(len_max) if len_max is not None else max_len
    dmin = int(dist_min) if dist_min is not None else min_dist
    dmax = int(dist_max) if dist_max is not None else max_dist
    if lmin > lmax: lmin, lmax = lmax, lmin
    if dmin > dmax: dmin, dmax = dmax, dmin

    df = df[df["ltr_len"].between(lmin, lmax, inclusive="both")]
    df = df[df["dist_to_tss"].between(dmin, dmax, inclusive="both")]

    # Sorting
    if sort_by:
        by = [s["column_id"] for s in sort_by]
        asc = [s["direction"] == "asc" for s in sort_by]
        df_sorted = df.copy()
        for c in by:
            if pd.api.types.is_categorical_dtype(df_sorted[c]):
                df_sorted[c] = df_sorted[c].astype(str)
        df = df_sorted.sort_values(by=by, ascending=asc, kind="mergesort")
    else:
        # Default: chromosomes first (chr1..chr22, chrX, chrY, chrM/MT) in natural order, then scaffolds; then by start
        df_sorted = df.copy()
        chrom = df_sorted["chrom"].astype(str)

        # mark true chromosomes
        is_chr = chrom.str.match(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")

        # numeric rank for chromosomes (1..22,X=23,Y=24,M/MT=25); NaN for scaffolds
        token = chrom.str.replace("chr", "", regex=False)
        num_raw = pd.to_numeric(token, errors="coerce")
        num_map = token.replace({"X": "23", "Y": "24", "M": "25", "MT": "25"})
        num_chr = pd.to_numeric(num_map, errors="coerce")

        # scaffold flag: chromosomes first (False=0), scaffolds after (True=1)
        df_sorted["_is_scaffold"] = ~is_chr
        # chromosome order (fill scaffolds with a large number so they sort after true chromosomes)
        df_sorted["_chr_ord"] = num_chr.fillna(1e9)

        df = df_sorted.sort_values(
            by=["_is_scaffold", "_chr_ord", "ltr_start"],
            ascending=[True, True, True],
            kind="mergesort",
        ).drop(columns=["_is_scaffold", "_chr_ord"])


    # Paging
    total = len(df)
    page_count = max(1, (total + page_size - 1) // page_size)
    start = page_current * page_size
    end   = start + page_size
    
    ltr_filtered = df.iloc[start:end].copy()
    ltr_filtered["Locus_display"] = ltr_filtered.apply(
        lambda r: fmt_coord(
            r.get("chrom"),
            r.get("ltr_start"),
            r.get("ltr_end"),
            r.get("ltr_strand")
        ),
        axis=1
    )

    return ltr_filtered.to_dict("records"), page_count

# From LTR table -> jump IGV
@app.callback(
    Output("genome-browser", "locus", allow_duplicate=True),
    Output("store-locus", "data", allow_duplicate=True),   # <-- add allow_duplicate
    Input("ltr-table", "selected_rows"),
    State("ltr-table", "data"),
    prevent_initial_call=True,
)
def jump_igv_ltr(selected_rows, table_data):
    from dash.exceptions import PreventUpdate
    import pandas as pd
    if not table_data or not selected_rows:
        raise PreventUpdate
    row = table_data[selected_rows[0]]
    chrom, start, end = row.get("chrom"), row.get("ltr_start"), row.get("ltr_end")
    if not chrom or pd.isna(start) or pd.isna(end):
        raise PreventUpdate
    pad = max(200, int(0.05 * (end - start)))
    locus = f"{chrom}:{max(1, int(start - pad))}-{int(end + pad)}"
    return locus, locus

@callback(
    Output("store-linked-internal-names", "data", allow_duplicate=True),
    Output("txt-link-summary", "children", allow_duplicate=True),
    Input("btn-link-ltrs", "n_clicks"),
    State("chk-link-use-selected", "value"),
    State("tbl-loci", "derived_virtual_data"),  # all rows after filter/sort
    State("tbl-loci", "data"),                  # current page (fallback)
    State("tbl-loci", "selected_rows"),
    prevent_initial_call=True,
)
def collect_internal_names(n_clicks, use_selected, dv_data, page_data, selected_rows):
    import pandas as pd
    # choose rows
    if use_selected and selected_rows:
        src = page_data or []
        rows = [src[i] for i in selected_rows if i is not None and i < len(src)]
        mode = "selected"
    else:
        rows = dv_data or page_data or []
        mode = "filtered"

    if not rows:
        return [], "No internal regions found from the current selection/filters."

    df = pd.DataFrame(rows)
    # The domains table uses "Locus" as the identifier column
    if "Locus" not in df.columns:
        return [], 'Column "Locus" not found in the domains table.'
    names = (
        df["Locus"].dropna().astype(str).drop_duplicates().sort_values().tolist()
    )
    return names, f"Linked {len(names)} internal region(s) from {mode} rows."


@callback(
    Output("store-linked-internal-names", "clear_data"),
    Output("txt-link-summary", "children", allow_duplicate=True),
    Output("tbl-loci", "selected_rows"),
    Input("btn-unlink", "n_clicks"),
    Input("btn-clear", "n_clicks"),
    prevent_initial_call=True,
)
def unlink_or_clear(n_unlink, n_clear):
    # Triggered by either Unlink or Clear Filters:
    if ctx.triggered_id in ("btn-unlink", "btn-clear"):
        # clear the link store, wipe the mini status text, and unselect rows in tbl-loci
        return True, "", []
    return no_update, no_update, no_update

# ---------- Clear filters ----------
@app.callback(
    Output("flt-subfamily", "value"),
    Output("flt-domain", "value"),
    Output("flt-mincov", "value"),
    Output("flt-ltr-subfamily", "value"),
    Output("flt-ltr-nmotifs", "value"),
    Output("flt-has-domain", "value"),
    Output("flt-ltrstatus", "value"),
    # NEW:
    Output("flt-u3r-class", "value"),
    Output("flt-u3r-feature", "value"),
    Output("flt-u3r-min-score", "value"),
    Output("flt-u3r-hide-lowconf", "value"),
    Input("btn-clear", "n_clicks"),
    prevent_initial_call=True,
)
def clear_filters(n):
    return (
        None,  # flt-subfamily
        None,  # flt-domain
        0.0,   # flt-mincov
        None,  # flt-ltr-subfamily
        0,     # flt-ltr-nmotifs
        False, # flt-has-domain
        None,  # flt-ltrstatus
        # U3R filters:
        None,  # flt-u3r-class
        None,  # flt-u3r-feature
        None,  # flt-u3r-min-score
        False, # flt-u3r-hide-lowconf
    )



# ---------- Link LTRs → U3/R/U5 table ----------
@callback(
    Output("store-linked-ltr-seqnames", "data", allow_duplicate=True),
    Output("txt-link-u3r-summary", "children", allow_duplicate=True),
    Input("btn-link-u3r", "n_clicks"),
    State("chk-link-u3r-use-selected", "value"),
    State("ltr-table", "derived_virtual_data"),  # filtered rows
    State("ltr-table", "data"),                  # current page (fallback)
    State("ltr-table", "selected_rows"),
    prevent_initial_call=True,
)
def collect_ltr_seqnames(n_clicks, use_selected, dv_data, page_data, selected_rows):
    import pandas as pd
    from dash.exceptions import PreventUpdate

    if not n_clicks:
        raise PreventUpdate

    if use_selected and selected_rows:
        src = page_data or []
        rows = [src[i] for i in selected_rows if i is not None and i < len(src)]
        mode = "selected"
    else:
        rows = dv_data or page_data or []
        mode = "filtered"

    if not rows:
        return [], "No LTRs found from the current selection/filters."

    df = pd.DataFrame(rows)
    if "sequence_name" not in df.columns:
        return [], 'Column "sequence_name" not found in the LTR table.'

    seqnames = (
        df["sequence_name"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    return seqnames, f"Linked {len(seqnames)} LTR(s) to U3/R/U5 table from {mode} rows."


@callback(
    Output("store-linked-ltr-seqnames", "clear_data"),
    Output("txt-link-u3r-summary", "children", allow_duplicate=True),
    Input("btn-unlink-u3r", "n_clicks"),
    prevent_initial_call=True,
)
def unlink_u3r_table(n_clicks):
    from dash.exceptions import PreventUpdate
    if not n_clicks:
        raise PreventUpdate
    return True, ""

@callback(
    Output("u3r-table", "data"),
    Input("store-linked-ltr-seqnames", "data"),
    Input("flt-u3r-class", "value"),
    Input("flt-u3r-feature", "value"),
    Input("flt-u3r-hide-lowconf", "value"),
    Input("flt-u3r-min-score", "value"),
)
def update_u3r_table(linked_seqnames, classes, feats, hide_lowconf, min_score):
    """
    Query the U3/R/U5 + PBS/PPT + signal annotations on demand via DuckDB.

    - If no LTRs are linked -> return empty table.
    - Otherwise, filter by linked LTR sequence_name and UI filters.
    - Hard-cap number of rows to avoid performance issues.
    """
    from math import isnan

    # Only show anything when there is an active link from the LTR table
    if not linked_seqnames:
        return []

    # Normalise inputs
    linked_seqnames = [str(s) for s in (linked_seqnames or [])]
    classes = classes or []
    feats = feats or []

    # Base WHERE clause: restrict to linked LTRs
    where_clauses = []
    params = [U3R_PARQUET]

    # sequence_name IN (?,?,...)
    placeholders_seq = ",".join(["?"] * len(linked_seqnames))
    where_clauses.append(f"sequence_name IN ({placeholders_seq})")
    params.extend(linked_seqnames)

    # Feature class filter (SEGMENT / PBS_PPT / SIGNAL)
    if classes:
        placeholders_cls = ",".join(["?"] * len(classes))
        where_clauses.append(f"feature_class IN ({placeholders_cls})")
        params.extend(classes)

    # Feature filter (U3, R, U5, PBS, PPT, TSS, Inr, PAS, ...)
    if feats:
        placeholders_feat = ",".join(["?"] * len(feats))
        where_clauses.append(f"feature IN ({placeholders_feat})")
        params.extend(feats)

    # Hide LOW_CONF segments (keep NULL / non-LOW_CONF)
    if hide_lowconf:
        where_clauses.append("(segment_conf IS NULL OR segment_conf <> 'LOW_CONF')")

    # Min score
    if min_score is not None and not (isinstance(min_score, float) and isnan(min_score)):
        try:
            thr = float(min_score)
            where_clauses.append("TRY_CAST(score AS DOUBLE) >= ?")
            params.append(thr)
        except Exception:
            pass

    # Assemble WHERE clause
    where_sql = " AND ".join(where_clauses) if where_clauses else "TRUE"

    # Limit rows sent to the browser
    MAX_ROWS = 1000

    sql = f"""
        SELECT
            feature,
            feature_class,
            feature_locus,
            score,
            segment_conf,
            sequence_name,
            subfamily,
            type,
            int_subfamily,
            nearest_gene,
            dist_to_tss
        FROM read_parquet(?)
        WHERE {where_sql}
        LIMIT {MAX_ROWS}
    """

    df = u3r_con.execute(sql, params).df()

    return df.to_dict("records")








# --------------------
# Run
# --------------------
if __name__ == "__main__":
    app.run(debug=True, dev_tools_hot_reload=False)
