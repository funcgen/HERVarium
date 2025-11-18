# Helper Scripts for HERVarium

This repository includes a small collection of preprocessing scripts used to generate the genomic tracks and precomputed tables displayed in the HERVarium web application.
These scripts are provided for transparency and reproducibility.
Most users will not need to run them, as all resulting bigBed files and precomputed tables are already distributed through the project’s Zenodo archive.

Below is a short description of each script.

## Python
```
prep_hervarium_tables.py
```
Generates all precomputed Parquet tables used by the HERVarium web-app. It merges LTR annotations, internal regions, domain calls, U3/R/U5 segments, PBS/PPT annotations, promoter/PAS signals, and TSS distances into unified tables stored under assets/precomputed/. This script is only required if rebuilding the full annotation tables from scratch.

## Shell scripts — building IGV/Dash tracks (bigBed)

These scripts convert BED files into bigBed tracks used by the IGV.js browser embedded in HERVarium.
They perform formatting, name simplification, score normalization, coordinate sorting, and the final bedToBigBed conversion.

``` 
make_ltr_bigbed.sh
``` 
- Cleans and converts the merged LTR BED file into a bigBed track (.bb).
- Ensures BED6 compliance (name, score, strand) and produces ERV_ltr_*.bb for visualization.

``` 
make_segments_bigbed.sh
``` 
- Processes U3/R/U5 LTR segments into a bigBed track.
- Standardizes segment names (U3, R, U5; HIGH/LOW_CONF), sorts, and builds the .bb file.

```
make_signals_bigbed.sh
```
- Generates a bigBed for promoter and PAS signals (TSS, Inr, DCE, PAS, cleavage motifs).
- Simplifies feature names and prepares HERV_LTR_U3_R_U5_signals.bb.

``` 
make_pbs_ppt_bigbed.sh
```
- Builds bigBed tracks for PBS and PPT elements detected near 5′ and 3′ LTRs.
- Cleans names (PBS/PPT), normalizes scores, sorts, and produces HERV_LTR_U3_R_U5_PBS_PPT.bb.

``` 
gencode_to_bigbed.sh
``` 
- Downloads GENCODE (user-specified release), converts the GTF into BED, replaces transcript IDs with gene symbols, sorts coordinates, and creates a GENCODE gene-symbol bigBed track for the genome browser.
- Used to generate the reference gene annotation track.

```
simplify_internal_and_ltr_names.sh
``` 
- Utility to clean and standardize BED name fields for internal regions and LTRs.
- Converts verbose RepeatMasker-style identifiers (e.g., prefix_pos_chr_start_end_strand_+) into compact prefix_chr_start_end_+ identifiers consistent with HERVarium.

## Usage

All scripts are optional and provided as reproducible components of the HERVarium annotation pipeline.
You do not need to run these scripts to use the web application—the corresponding bigBed tracks and precomputed tables can be downloaded directly from the HERVarium Zenodo repository.
