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
- `signature_dir = pl.Path("/add/path/here/signatures_canceronly/")`: path to where the cancer-specific signatures computed in the analysis notebook [`6. cancer-specific-signature.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/6.%20cancer-specific-signature.ipynb) were saved.
- `gencode = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `bulk = pd.read_csv("/add/path/here/LUD2015-005_RNAseq_featureCounts.tsv",sep="\t",index_col=0)`: path to where the Carroll et al. bulk dataset was saved. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)

## 4. `compare-Nowicki-BE.ipynb`

- `signature_dir = pl.Path("/add/path/here/")`: path to where the cNMF program marker genes were saved in the analysis notebook [`5. cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/5.%20cNMFCancerCells-perPatient.ipynb).
- `science_path = pl.Path("/add/path/here/")`: path to where the signatures derived from the Nowicki-Osuch et al. paper are located. (**DWNL**)

## 5. `runpySCENIC-external part 1.ipynb`

- `adata = sc.read_h5ad("/add/path/here/Carroll_EAC_raw.h5ad")`: path to where the Carroll et al. single-cell dataset was saved. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)
- `f_loom_path_scenic = "/add/path/here/adata_filtered_scenic.loom"`: path to where the Carroll et al. counts will be saved as a loom file - input to the pySCENIC pipeline.
- `adata_croft = sc.read_h5ad("/add/path/here/GSE222078_adata.h5ad")`: path to where the processed data from Croft et al. was saved in [`Croft-validation-set.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/2.%20Croft-validation-set.ipynb). (**DWNL**)
- `f_loom_path_scenic = "/add/path/here/adata_croft_filtered_scenic.loom"`: path to where the Croft et al. counts will be saved as a loom file - input to the pySCENIC pipeline.
- `adata_luo = sc.read_h5ad("/add/path/here/GSE210347_fibroblast_counts.h5ad")`: path to where the processed data from Luo et al. was saved in [`CAFatlas-fibroblast-validation.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/1.%20CAFatlas-fibroblast-validation.ipynb). (**DWNL**)
- `f_loom_path_scenic = "/add/path/here/adata_luo_filtered_scenic.loom"`: path to where the Luo et al. counts will be saved as a loom file - input to the pySCENIC pipeline. 

#### RUN PYSCENIC

The SCENIC method needs to be run on the three external datasets at this point. Further instructions can be found [here](link).

## 5. `runpySCENIC-external part 2.ipynb`

- `orig_corrs[state] = pd.read_csv(f"/add/path/here/{state}_score_triad_corr.csv",index_col=0)`: path to where the correlation results were saved in the analysis notebook [`9. SCENICplus_analyze_cNMF.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/9.%20SCENICplus-analyze-cNMF.ipynb).
- `adata = sc.read_h5ad("/add/path/here/Carroll_EAC_raw.h5ad")`: path to where the Carroll et al. single-cell dataset was saved. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)
- `lf = lp.connect("/add/path/here/pyscenic_carroll_output.loom", mode='r+', validate=False )`: path to where the results from the SCENIC analysis on the Carroll et al. dataset were saved.
- `signature_dir = pl.Path("/add/path/here")`: path to where the cNMF program marker genes were saved in the analysis notebook [`5. cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/5.%20cNMFCancerCells-perPatient.ipynb).
- `toptfs = pd.read_csv("/add/path/here/toptfs_top20.csv",index_col=0)`: path to where the top 20 candidate mTFs were saved in the analysis notebook [`SCENICplus-analyze-cNMF.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/9.%20SCENICplus-analyze-cNMF.ipynb)
- `all_tfs = pd.read_csv("/add/path/here/utoronto_human_tfs_v_1.01.txt",header=None).values.ravel()`:  path to the list of TFs provided as an auxiliary SCENIC file (see SCENIC/SCENIC+ instructions for more info). (**DWNL**)
- `adata_croft = sc.read_h5ad("/add/path/here/GSE222078_adata.h5ad")`: path to where the processed data from Croft et al. was saved in [`Croft-validation-set.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/2.%20Croft-validation-set.ipynb). (**DWNL**)
- `celltypes = pd.read_csv("/add/path/here/highLevelCellTypes.csv",index_col=0)`: path to where the cell types of the Croft data are located.
- `lf = lp.connect("/add/path/here/pyscenic_croft_output.loom", mode='r+', validate=False )`: path to where the results from the SCENIC analysis on the Croft et al. dataset were saved.
- `adata_luo = sc.read_h5ad("/add/path/here/GSE210347_fibroblast_counts.h5ad")`: path to where the processed data from Luo et al. was saved in [`CAFatlas-fibroblast-validation.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/1.%20CAFatlas-fibroblast-validation.ipynb). (**DWNL**)
- `lf = lp.connect("pyscenic-external-results/pyscenic_luo_output.loom", mode='r+', validate=False )`: path to where the results from the SCENIC analysis on the Luo et al. dataset were saved.
- `signature_dir2 = pl.Path("/add/path/here/")`: path to where the fibroblast marker genes were saved in the analysis notebook [`3.plotting-TME-subsets.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/3.%20plotting-TME-subsets.ipynb)

## 6. `Hoefnagel-bulk-validation.ipynb`

- `gex_df = pd.read_csv("/add/path/here/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt",sep="\t").iloc[1:,:].T`: path to where the Hoefnagel et al. was downloaded (**DWNL**).
- `gencode_mapping = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `survival_clin = pd.read_csv("/add/path/here/data.SPSS.subselect.txt", sep="\t", index_col=0)`: path to where the file with the survival information was saved (**DWNL**).
- `gse = GEOparse.get_GEO(geo="GSE207526", destdir="/add/path/here/")`: path to where the clinical information accessible through GEOparse will be saved.
- `signature_dir = pl.Path("/add/path/here/signatures_canceronly/")`: path to where the cancer-specific signatures computed in the analysis notebook [`6. cancer-specific-signature.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/6.%20cancer-specific-signature.ipynb) were saved.

## 7. `TCGA-bulk-validation.ipynb`

- `gencode_mapping = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `eac_tcga_dir = pl.Path("/add/path/here/")`: path to the directory where all the TCGA-related files linked in the main [README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md) are located.
- `signature_dir = pl.Path("/add/path/here/signatures_canceronly")`: path to where the cancer-specific signatures computed in the analysis notebook [`6. cancer-specific-signature.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/6.%20cancer-specific-signature.ipynb) were saved.

#### NEED TO RUN BAYESPRISM

At this point, you need to have run the R scripts of BayesPrism (see more instructions [here](https://github.com/vanallenlab/EAC-multiome/blob/main/code/R/scripts/BayesPrism/README.md). 

##8. `Link-States-Ecotypes-BayesPrism.ipynb`

- `tcga_dir = pl.Path("/add/path/here/")`: path to the directory where all the TCGA-related files linked in the main [README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md) are located.
Same as for the [`7. TCGA-bulk-validation.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/7.%20TCGA-bulk-validation.ipynb) script.
- `bp_eac = pd.read_csv("/add/path/here/eac_purity.csv",index_col=0)`: this is the path to the results of BayesPrism for TCGA deconvolution.
- `bp_eac_gse = pd.read_csv("/add/path/here/eac_gse_purity.csv",index_col=0)`: this is the path to the results of BayesPrism for Hoefnagel et al. deconvolution
- `bp_eac_carroll = pd.read_csv("/add/path/here/eac_carroll_purity.csv",index_col=0)`: this is the path to the results of BayesPrism for Carroll et al. deconvolution
- `purity = pd.read_csv("/add/path/here/TCGA_absolute_purity.txt",index_col=0,sep="\t")` this is the path to the ABSOLUTE estimated purity of TCGA samples. (**DWNL**)
- `gex_df = pd.read_csv("/add/path/here/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt",sep="\t").iloc[1:,:].T`: path to where the Hoefnagel et al. was downloaded (**DWNL**).
- `gencode_mapping = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`: this file is a subset of the Gencode annotation v41 that contains "seqname", "gene_name", "start", "end", "strand", and "gene_id". (**DWNL**)
- `gex_df2 = pd.read_csv("/add/path/here/bulk_preprocessed.csv",index_col=0).T`: this is the path to where the Carroll et al. bulk dataset **after preprocessing** was saved. In practice, this means the LUD2015-005_RNAseq_featureCounts.tsv was simply transformed to keep only the genes on the rows and the patients on the columns. *Note: this data may vary depending on the exact processing performed on the Carroll et al. dataset* (**DWNL**)
- `signature_dir = pl.Path("/add/path/here/signatures_canceronly/")`: path to where the cancer-specific signatures computed in the analysis notebook [`6. cancer-specific-signature.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/6.%20cancer-specific-signature.ipynb) were saved.
