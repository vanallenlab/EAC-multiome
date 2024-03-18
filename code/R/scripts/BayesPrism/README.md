# EAC-multiome

This folder contains all scripts used to run BayesPrism to deconvolve the bulk data from the TCGA, Hoefnagel et al. and Carroll et al. cohorts.

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be download prior to running the scrips is indicated by the <span style="color:red">DWNL</span> tag.

## `run_BPrism.R`

- `geneMap <- "/add/path/here"`: this is the path to the gene mapping to map ENCODE to official gene ID <span style="color:red">DWNL</span>
- `path <- "add/path/here"`: path to the general folder where the data is located
- `seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))`: path to single-cell Carroll et al. data <span style="color:red">DWNL</span>
- `data <- read.table(paste0(path, expr), sep=',', row.names=1, header=TRUE)`: path to the bulk TCGA raw counts <span style="color:red">DWNL</span>
- `save(bp.res, file=paste0(path, "/add/path/here/eac.bp.res.RData"))`: path to where the results should be saved

The information can easily be extracted from the `.RData` object using the following commands
```
# to obtain the proportions of each cell type
theta <- get.fraction (bp=bp.res, which.theta="final",
             state.or.type="type")
```
 