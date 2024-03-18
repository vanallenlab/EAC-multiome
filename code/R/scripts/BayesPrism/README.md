# EAC-multiome

This folder contains all scripts used to run BayesPrism to deconvolve the bulk data from the TCGA, Hoefnagel et al. and Carroll et al. cohorts.

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

## `run_BPrism.R`

- `geneMap <- "/add/path/here"`: this is the path to the gene mapping to map ENCODE to official gene ID (**DWNL**)
- `path <- "add/path/here"`: path to the general folder where the data is located
- `seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))`: path to single-cell Carroll et al. data (**DWNL**)
- `expr <- "add/path/here/TCGA-ESCA.htseq_counts.tsv"`: path to the bulk TCGA raw counts (**DWNL**)
- `save(bp.res, file=paste0(path, "/add/path/here/eac.bp.res.RData"))`: path to where the results should be saved

## `run_BPrism_Hoefnagel.R`

- `path <- "add/path/here"`: path to the general folder where the data is located
- `seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))`: path to single-cell Carroll et al. data (**DWNL**)
- `expr <- "/add/path/here/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt"`: path to the downloaded Hoefnagel et al raw counts (**DWNL**). 
- `save(bp.res, file=paste0(path, "/add/path/here/eac-gse207526.bp.res.RData"))`: path to where the results should be saved

  ## `run_BPrism_Carroll.R`

- `path <- "add/path/here"`: path to the general folder where the data is located
- `seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))`: path to single-cell Carroll et al. data (**DWNL**)
- `expr <- "add/path/here/bulk_preprocessed.csv"`: path to the downloaded Carroll et al raw counts (**DWNL**). *Note: This file has been preprocessed to have official GENE ID on the rows and patient names on the columns.*
- `save(bp.res, file=paste0(path, "/add/path/here/eac-carroll.bp.res.RData"))`: path to where the results should be saved

The information can easily be extracted from the `.RData` object using the following commands
```
# to obtain the proportions of each cell type
theta <- get.fraction (bp=bp.res, which.theta="final",
             state.or.type="type")
``` 
