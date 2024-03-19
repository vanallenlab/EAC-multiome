# EAC-multiome

This folder contains all the notebooks necessary to preprocess the snRNA data on a per-sample basis. 

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

## Step 0

The input files to this part are already filtered for ambient RNA and doublets using Cellbender and Scrublet, with default parameters for both methods. If you wish to reprocess the data from the raw files, you will need to run [Cellbender](https://github.com/broadinstitute/CellBender) and [Scrublet](https://github.com/swolock/scrublet) first. Otherwise, processed .h5ad are available for download [here](link). #TBD

All the files require the same input, described here
## `Aguirre_XXX.ipynb` or `CCG1153_XXX.ipynb`

- `datapath = pl.Path("/add/path/here/")`: path to the directory where the Cellbender and Scrublet filtered sample-level data are located. (**DWNL**)
- `gencode_df = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0).set_index("gene_name")`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `infercnv_annot_dir = pl.Path("/add/path/here/infercnv-annotations")`: this is the directory where the inferCNV annotations are saved.
- `resdir = pl.Path("/add/path/here/")`: this is the directory where the preprocessed sample-level files are saved.


