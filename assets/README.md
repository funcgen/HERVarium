# HERVarium – Assets Directory

This directory contains the **static and data files** required by the HERVarium web-application.  

Only a small subset of files (logos) is tracked in the GitHub repository. **All functional data files used by the web-app must be downloaded from the Zenodo releases** of the HERVarium project.

## Contents

### 1. Logos and css (included in GitHub)
These files are the only assets shipped directly with the repository:

- `logo_cnag.jpg`
- `logo_generalitat.png`
- `logo_eu.png`
- `styles.css`

### 2. BigBed Genome Browser Tracks (download from Zenodo)
HERVarium uses several `.bb` tracks for visualization in the embedded IGV.js browser:

- Internal HERV regions  
- Conserved domain calls  
- Merged LTRs  
- U3/R/U5 segments  
- PBS and PPT annotations  
- Promoter and PAS signal positions  
- FIMO TFBM hits  
- GENCODE gene symbols  

These **are not included in the GitHub repository**, due to their size.  
Download the complete collection from Zenodo:

👉 **Zenodo: [Domain-Level Annotations and Conservation Scores for Human Endogenous Retroviruses](https://doi.org/10.5281/zenodo.16318927)**  
👉 **Zenodo: [Structural Annotation and Transcription Factor Motif Maps for Human Endogenous Retroviruses](https://doi.org/10.5281/zenodo.17602210)**  


Place all downloaded `.bb` files inside this `assets/` directory.

### 3. Reference Genome Files (download from Zenodo)
HERVarium requires the GRCh38 primary assembly FASTA and its index:

- `GRCh38.primary_assembly.genome.fa`
- `GRCh38.primary_assembly.genome.fa.fai`

These files are also hosted on Zenodo and must be placed in this directory.

### 4. Precomputed Annotation Tables
The `assets/precomputed/` directory contains the Parquet and JSON metadata files used by the web-app for instant filtering.  
These files are **not stored in GitHub**.  
See the `assets/precomputed/README.md` for detailed instructions.

## After Downloading the Assets
After placing the required files in this directory, HERVarium can be run locally via:

```
python app.py
```

The web-app will automatically detect the bigBed tracks and precomputed tables.
