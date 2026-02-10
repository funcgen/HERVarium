# HERVarium  
### Genome-wide atlas of HERV internal domains, LTR regulatory motifs, and U3–R–U5 architecture

HERVarium is an interactive Dash-based web application to explore internal protein-coding
domains and LTR regulatory features of human endogenous retroviruses (HERVs).  
It integrates two genome-wide resources:

1. **Internal retroviral protein domains (GyDB/DFAM HMM profiles)**  
2. **LTR regulatory architecture**, including:
   - U3–R–U5 segments  
   - Promoter motifs (TATA, Inr, DPE, MTE, DCE…)  
   - PAS hexamers  
   - PBS (tRNA binding sites)  
   - PPT (polypurine tracts)  
   - ∼900 transcription factor binding motifs (TFBMs; FIMO)  

HERVarium allows you to browse loci through an embedded **IGV browser**, apply biological filters,
and link internal regions → LTRs → U3/R/U5 regulatory features.

## 🧬 Data availability

All raw and derived datasets used in HERVarium are publicly available in two Zenodo repositories:

- **LTR regulatory atlas (U3–R–U5, motifs, PBS/PPT)**  
  https://doi.org/10.5281/zenodo.17602210

- **Internal domain annotation (HERV ORFs + GyDB HMMs)**  
  https://doi.org/10.5281/zenodo.16318927

## 📁 Repository structure

```
HERVarium/
│
├── app.py                          # Main Dash application
├── environment.yml                 # Reproducible conda environment
│
├── scripts/                        # Helper scripts to rebuild HERVarium assets
│   ├── prep_hervarium_tables.py
│   ├── gencode_to_bigbed.sh
│   ├── make_ltr_bigbed.sh
│   ├── make_segments_bigbed.sh
│   ├── make_signals_bigbed.sh
│   ├── make_pbs_ppt_bigbed.sh
│   ├── convert_fimo_to_bigbed.sh
│   ├── simplify_fimo_bed_name.py
│   ├── simplify_domains_bed.sh
│   ├── simplify_internal_and_ltr_names.sh
│   ├── make_gtex_expressed.sh
│   ├── build_motif_duckdb.py
│   └── change_bb_colors.sh
│
├── bin/
│   └── bedToBigBed                 # UCSC bedToBigBed binary (vendored copy)
│
├── assets/
│   ├── genome/
│   │   ├── GRCh38.primary_assembly.genome.fa
│   │   ├── GRCh38.primary_assembly.genome.fa.fai
│   │   └── GRCh38.primary_assembly.genome.tar.xz
│   │
│   ├── gencode/
│   │   ├── gencode.v48.primary_assembly.annotation.gtf
│   │   ├── gencode.v48.bed
│   │   ├── gencode.v48.genepred
│   │   ├── gencode.v48.genesymbols.bed
│   │   ├── gencode.v48.genesymbols.sorted.bed
│   │   └── gencode.v48.genesymbols.bb
│   │
│   ├── internals/
│   │   ├── ERV_full_plus_components.bed
│   │   ├── ERV_full_plus_components.map.tsv
│   │   ├── ERV_GyDB_v6_domains.bed
│   │   ├── HERV_internal_v6.bed
│   │   ├── HERV_internal_simplified.bed
│   │   ├── HERV_internal_domains_simplified.bed
│   │   ├── HERV_loci_annotated_domains.tsv
│   │   └── INTERNAL_fully_annotated.tsv
│   │
│   ├── ltr/
│   │   ├── ERVs_LTRs_merged_v4.bed
│   │   ├── ERVs_LTRs_merged_v4.simplified.bed
│   │   ├── ERV_ltr_v4_merged.simplified.clean.bed
│   │   ├── ERV_ltr_v4_merged.simplified.sorted.bed
│   │   ├── ERV_ltr_merged.simplified.bb
│   │   ├── LTR_fully_annotated.tsv
│   │   │
│   │   ├── segments/
│   │   │   ├── HERV_LTR_U3_R_U5_catalogue.tsv
│   │   │   ├── HERV_LTR_U3_R_U5_segments_allconf.bed
│   │   │   ├── HERV_LTR_U3_R_U5_segments_allconf.clean.bed
│   │   │   ├── HERV_LTR_U3_R_U5_segments_allconf.sorted.bed
│   │   │   ├── HERV_LTR_U3_R_U5_segments_allconf.bb
│   │   │   ├── HERV_U3_R_U5_segments_highconf.bed
│   │   │   ├── HERV_LTR_U3_R_U5_signals.bed
│   │   │   ├── HERV_LTR_U3_R_U5_signals.clean.bed
│   │   │   ├── HERV_LTR_U3_R_U5_signals.sorted.bed
│   │   │   └── HERV_LTR_U3_R_U5_signals.bb
│   │   │
│   │   └── tfbm/
│   │       ├── fimo_parsed_v4.tsv
│   │       ├── fimo_parsed_v4.bed
│   │       ├── fimo_parsed_v4.sorted.bed
│   │       ├── fimo_parsed_v4.fixed.bed
│   │       └── fimo_parsed.bb
│   │
│   ├── precomputed/
│   │   ├── agg.parquet
│   │   ├── ltr.parquet
│   │   ├── ltr_u3r_u5.parquet
│   │   ├── domains_meta.json
│   │   ├── ltr_meta.json
│   │   └── ltr_u3r_u5_meta.json
│   │
│   ├── hg38.chrom.sizes
│   ├── styles.css
│   ├── favicon.ico
│   └── logos/
│       ├── hervarium_logo.png
│       ├── logo_cnag.jpg
│       ├── logo_generalitat.png
│       └── logo_eu.png
│
└── README.md

```

Each Zenodo record contains:  
• BED/BigBed files  
• FASTA files  
• Tables (TSV/Parquet)  
• Metadata JSON  
• Documentation of file formats  


## ⚙️ Installation

HERVarium can be installed and run locally using **conda**.

### 1. Clone the repository

```
git clone https://github.com/<YOUR-USERNAME>/HERVarium.git
cd HERVarium
``` 

### 2. Create the environment

#### 📦 environment.yml

Below is a complete conda environment tested with HERVarium:

```
name: hervarium
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.10
  - pip
  - dash
  - dash-bio
  - dash-bootstrap-components
  - flask-caching
  - pandas
  - numpy
  - duckdb
  - pyarrow
  - gunicorn        # optional for deployment
  - aiohttp         # optional igv.js compatibility
  - pip:
      - dash==2.14.2
      - dash-bio==1.0.2
      - dash-bootstrap-components==1.6.0
      - flask-caching
``` 


``` 
conda env create -f environment.yml
conda activate hervarium
```

### 3. Download the precomputed annotation files

Download the Zenodo HERVarium data bundle:
- [https://doi.org/10.5281/zenodo.18551737](https://doi.org/10.5281/zenodo.18551737)

From each file, copy them to their corresponding folder. Keep the file names unchanged (the app expects these exact names). 

## ▶️ Run HERVarium locally

```
python app.py
```

Then open:
```
http://127.0.0.1:8050
```
## 🧬 Usage

Main functionalities:

### 1. Locus Browser (IGV)

- Navigate genomic coordinates
- View GENCODE annotations
- View HERV internal regions, LTRs, U3/R/U5 segments, motifs, and TFBMs
- Optional ENCODE DNase tracks (cell-type selectable)
- Optional GTEx RNA-seq tracks (tissue selectable)

### 2. Internal domain table

- Filter by subfamily, domain class, coverage, LTR status
- Link internal regions → corresponding LTRs
- Export results to CSV

### 3. LTR regulatory table

- Filter by subfamily, LTR type, distance to TSS, #motifs
- Link LTRs → U3/R/U5 regulatory features

### 4. U3/R/U5 + PBS + PPT + signals

- Query dynamically via DuckDB
- Filter by feature class, feature name, min score, confidence
- Export results to CSV

### 📖 Citation

If you use HERVarium in your work, please cite:

>Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). Regulatory Features and Functional Specialization of Human Endogenous Retroviral LTRs: A Genome-Wide Annotation and Analysis via HERVarium. Manuscript in preparation.

Data citations:
> Internal domain annotation: https://doi.org/10.5281/zenodo.16318927

> LTR regulatory atlas: https://doi.org/10.5281/zenodo.17602210  

