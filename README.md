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

## 🚀 Quick start 

For reviewers and users who want to run HERVarium locally with minimal setup,
we provide a **self-contained data bundle** hosted on Zenodo.

This bundle includes:
- The complete HERVarium application
- All precomputed annotation assets (BED/BigBed, Parquet, FASTA)
- The conda environment file
- Directory structure expected by the app

**No manual data assembly is required.**

### Download and run HERVarium locally

1. Download the HERVarium data bundle:
   - https://doi.org/10.5281/zenodo.18551737

2. Unpack the archive:
   ```bash
   tar -xvf hervarium.tar.xz
   cd HERVarium
   ```
3. Create and activate the conda environment:

   ```bash
   conda env create -f environment.yml
   conda activate hervarium
   ```
4. Run the application:

   ```bash
   python app.py

5. Open your browser at:

   ```
   http://127.0.0.1:8050
   ``` 
This is the recommended installation method during the preprint and initial release phase.


## 🧬 Data availability

All raw and derived datasets used in HERVarium are publicly available in three Zenodo repositories:

- **LTR regulatory atlas (U3–R–U5, motifs, PBS/PPT)**  
  https://doi.org/10.5281/zenodo.17602210

- **Internal domain annotation (HERV ORFs + GyDB HMMs)**  
  https://doi.org/10.5281/zenodo.16318927

The recommended way to obtain a fully functional local installation of HERVarium
is via the prepackaged **HERVarium data bundle** hosted on Zenodo: https://doi.org/10.5281/zenodo.18551737.  
This archive contains all scripts, assets, and directory structure required to run
the application locally without additional downloads.
   

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
│   └── make_gtex_expressed.sh
│
├── bin/
│   └── bedToBigBed                 # UCSC bedToBigBed binary (vendored copy)
│
├── assets/
│   ├── genome/
│   │   ├── GRCh38.primary_assembly.genome.fa
│   │   └── GRCh38.primary_assembly.genome.fa.fai
│   │
│   ├── gencode/
│   │   └── gencode.v48.genesymbols.bb
│   │
│   ├── internals/
│   │   ├── HERV_internal_simplified.bed
│   │   ├── HERV_internal_domains_simplified.bed
│   │
│   ├── ltr/
│   │   ├── ERV_ltr_merged.simplified.bb
│   │   │
│   │   ├── segments/
│   │   │   ├── HERV_LTR_U3_R_U5_segments_allconf.bb
|   |   |   ├── HERV_LTR_U3_R_U5_PBS_PPT.bb
│   │   │   └── HERV_LTR_U3_R_U5_signals.bb
│   │   │
│   │   └── tfbm/
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

>Regulatory Features and Functional Specialization of Human Endogenous Retroviral LTRs: A Genome-Wide Annotation and Analysis via HERVarium. Tomàs Montserrat-Ayuso, Aurora Pujol, Anna Esteve-Codina. bioRxiv 2026.02.17.706328; doi: https://doi.org/10.64898/2026.02.17.706328

Data citations:
> Internal domain annotation: https://doi.org/10.5281/zenodo.16318927

> LTR regulatory atlas: https://doi.org/10.5281/zenodo.17602210  

