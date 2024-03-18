# EAC-multiome

This folder contains the scripts to preprocess the ATAC data. Data that should be produced through an intermediate script is indicated by (**DWNL, script-name**). 

## `signac-MACS2-callpeaks.R`
- `basedir <- "/add/path/here"`: path to where the per-sample raw ATAC data is present. This would be the output of the CellRanger ARC pipeline. We will use the `per_barcode_metrics.csv`, `filtered_feature_bc_matrix.h5`, and `atac_fragments.tsv.gz` files.
- `macsdir <- "/add/path/here"`: path to where the MACS2 output will be saved
- `macs2.path="/add/path/here/macs2"`: path to where the MACS2 executable is located (Instructions to installation [here](https://github.com/macs3-project/MACS/wiki/Install-macs2))

## `signac-MACS2-combine.R`
- `ct.annotations <- read.csv("add/path/here/adata_cNMF_scores_wtop.csv")` (**DWNL, script-name**)
- `basedir <- "/add/path/here/"`: path to where the per-sample raw ATAC data is present. This would be the output of the CellRanger ARC pipeline. We will use the `per_barcode_metrics.csv` and `atac_fragments.tsv.gz` files.
- `macsdir <- "/add/path/here/"`: path to where the MACS2 output from `signac-MACS2-callpeaks.R` was saved
- `respath <- "/add/path/here/"`: path to where the data should be saved 
