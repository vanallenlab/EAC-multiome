# EAC-multiome

This folder contains all notebooks used for the validation analyses using external cohorts. 

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

These notebooks should be run in the following order:

## 1. `CAFatlas-fibroblast-validation.ipynb`

- `signature_dir = pl.Path("/add/path/here/fibroblast/")`: path to where the fibroblast marker genes were saved in the analysis notebook [`3.plotting-TME-subsets.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/3.%20plotting-TME-subsets.ipynb)
- `metadata = pd.read_csv("/add/path/here/GSE210347_meta.txt.gz", sep="\t",index_col=0)`: path to where the Luo et al. pan-cancer CAF atlas metadata was saved. (**DWNL**)
- `adata = sc.read_h5ad("/add/path/here/GSE210347_counts.h5ad")`: path to where the where the Luo et al. pan-cancer CAF atlas count data was saved. (**DWNL**)
- `fibadata.write_h5ad("/add/path/here/GSE210347_fibroblast_counts.h5ad")`: path to where the data will be saved.

## 2. `Croft-validation-set.ipynb`

- `datadir = pl.Path("/add/path/here")`: path to where the original data downloaded from GEO is located.  (**DWNL**)
- `adata.write("/add/path/here/GSE222078_adata.h5ad")` and  `adata = sc.read_h5ad("/add/path/here/GSE222078_adata.h5ad")`: path to where the processed data will be saved.
- `signature_dir = pl.Path("/add/path/here")`: path to where the cNMF program marker genes were saved in the analysis notebook [`5. cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/5.%20cNMFCancerCells-perPatient.ipynb).
- `toptfs = pd.read_csv("/add/path/here/toptfs_top20.csv",index_col=0)`: path to where the top 20 candidate mTFs were saved in the analysis notebook [`SCENICplus-analyze-cNMF.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/9.%20SCENICplus-analyze-cNMF.ipynb)

## 3. `Carroll-validation-set.ipynb`

- `cell_cycle_genes = [x.strip() for x in open('/add/path/here/regev_lab_cell_cycle_genes.txt')]`: path to where the cell cycle genes are saved. (**DWNL**)
- `adata = sc.read_h5ad("/add/path/here/Carroll_EAC_raw.h5ad")`: path to where the Carroll et al. single-cell dataset was saved. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)
- `clinical = pd.read_csv("/add/path/here/carroll_clinical.csv", index_col=0)`: path to where the Carroll et al. clinical information was saved. (**DWNL**)
- `signature_dir = pl.Path("/add/path/here")`: path to where the cNMF program marker genes were saved in the analysis notebook [`5. cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/5.%20cNMFCancerCells-perPatient.ipynb).
- `signature_dir2 = pl.Path("/add/path/here/")`: path to where the fibroblast marker genes were saved in the analysis notebook [`3.plotting-TME-subsets.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/3.%20plotting-TME-subsets.ipynb)
- `toptfs = pd.read_csv("/add/path/here/toptfs_top20.csv",index_col=0)`: path to where the top 20 candidate mTFs were saved in the analysis notebook [`SCENICplus-analyze-cNMF.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/9.%20SCENICplus-analyze-cNMF.ipynb)
- `gencode = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `bulk = pd.read_csv("/add/path/here/LUD2015-005_RNAseq_featureCounts.tsv",sep="\t",index_col=0)`: path to where the Carroll et al. bulk dataset was saved. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)

## `external_GSE207527.ipynb`

- `gex_df = pd.read_csv("/add/path/here/GSE207526/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt",sep="\t").iloc[1:,:].T`
- `gencode_mapping = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`
- `gse = GEOparse.get_GEO(geo="GSE207526", destdir="/add/path/here")`
- `survival_clin = pd.read_csv("/add/path/here/data.SPSS.subselect.txt", sep="\t", index_col=0).set_index("FileName.GenomeScan")`
- `signature_dir = pl.Path("/add_path_here/")`



