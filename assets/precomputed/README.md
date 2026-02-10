# HERVarium – Precomputed Annotation Tables

This directory contains the **precomputed Parquet tables and JSON metadata** used by the HERVarium web-application. These files support fast filtering, table paging, and on-demand queries of LTRs, internal regions, and U3/R/U5 features. **No precomputed files are included in the GitHub repository.** They must be downloaded from the Zenodo releases associated with the project.

## File Overview

The following files must be placed in this directory:

### Parquet Tables
- `ltr.parquet`  
  Full annotation for all LTR elements (classification, distances, motif counts).

- `agg.parquet`  
  Aggregated summaries linking LTRs to internal regions, structural class, and domain conservation.

- `ltr_u3r_u5.parquet`  
  U3/R/U5 segments, PBS/PPT attributes, promoter/PAS signals, and confidence scores.

### Metadata Files (JSON)
- `ltr_meta.json`
- `domains_meta.json`
- `ltr_u3r_u5_meta.json`

These files specify valid categories, filter ranges, and table schemas used by the application.

## Downloading the Data

All files are available from the Zenodo archive:

👉 [https://doi.org/10.5281/zenodo.18551737](https://doi.org/10.5281/zenodo.18551737)  

Download the `.parquet` and `.json` files and place them inside this directory.


## Recreating the Tables (optional)

Advanced users can regenerate these tables using:

scripts/prep_hervarium_tables.py

This requires the full collection of annotation inputs (domains, LTRs, motifs, U3/R/U5 BEDs, etc.).

For most users, the **Zenodo-distributed tables are recommended**.

## Usage

When running:
```
python app.py
```

HERVarium automatically loads the Parquet tables from this directory to populate its interactive tables, filters, and genome browser selection logic.

---

