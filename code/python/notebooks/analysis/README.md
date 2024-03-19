# EAC-multiome

This folder contains all notebooks used for the analysis of the discovery cohort. 

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

These notebooks should be run in the following order:

## 1. `cohort-plot.ipynb`

- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**)
- `toy_comut.figure.savefig("add/path/here/full_cohort_plot.svg", dpi=200)`: path to where to save the Figure.

## 2. `tme-cleaning-analysis.ipynb`

- `datadir = pl.Path("add/path/here/celltype_h5ad_files/")`: path to the folder that contains the output of the preprocessing notebooks. Should contain one `.h5ad` file per patient.
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**)
- `refined_annotations.to_csv("add/path/here/refined_annotations.csv")`: path to where to save the refined annotations derived from this notebook. 
- `adata.write("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object will be saved. 

## 3. `plotting-TME-subsets.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**).

## 4. `run-milopy.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv",index_col=0)`: path to where the refined annotations derived from [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb) were saved. 

## 5. `cNMFCancerCells-perPatient.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv", index_col=0)`: path to where the refined annotations derived from [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb) are saved.
- `pd.Series(adata.var_names).to_csv("/add/path/here/eac_gene_names.csv")`: path to where to save the genes present in the dataset
- `patadata.write(f"/add/path/here/{sample}_subadata_cNMF.h5ad")` and `sample_file = f"/add/path/here/{sample}_subadata_cNMF.h5ad"`: path to where to save patient-level h5ad with all annotations - these are used to then derive the cNMF programs on a sample level. this must be replaced for each sample paragraph
- `marker_genes[cl].to_csv(f"/add/path/here/cNMF_{cl}.csv")`: path to where to save the marker genes for each cNMF program
- `("cNMF_"+ cluster_assignment.astype(str)).to_csv("/add/path/here/cNMF_program_assignment_cluster.csv")`: path to where to save the association between the derived programs and the cluster assignments
- `itay_MPs = pd.read_csv("/add/path/here/ItayTiroshHeterogeneityMPs.csv")`: path to where the gene programs derived from the Barkley et al. paper are located. (**DWNL**)
- `gene_sets='/add/path/here/h.all.v7.4.symbols.gmt',`:  path to where the .gmt for the MSigDB cancer hallmarks is located. (**DWNL**)

## 6. `cancer-specific-signature.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
-`marker_genes[cl] = pd.read_csv(f"/add/path/here/cNMF_{cl}.csv",index_col=0)`: path to where the marker genes were saved in [`cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cNMFCancerCells-perPatient.ipynb.ipynb).
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv",index_col=0)`: path to where the refined annotations derived from [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb) were saved. 
- `resdir = pl.Path("/add/path/here/")`: path to where to save the results

#### RUN THE SNATAC ANALYSIS HERE
The ATAC analysis R scripts need to be run here. See more info [here](https://github.com/vanallenlab/EAC-multiome/tree/main/code/R/scripts/preprocessing-snATAC).

## 7. `cohort-viz.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv", index_col=0)`: path to where the refined annotations derived from [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb) are saved.
- `adata.obs[["refined_annotations","sample_id"]].to_csv("/add/path/here/refined_annotations_wsampleid.csv")`: path to where the refined annotations+sample ID will be saved.
- `marker_genes[cl] = pd.read_csv(f"/add/path/here/cNMF_{cl}.csv",index_col=0)`: path to where the marker genes were saved in [`cNMFCancerCells-perPatient.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cNMFCancerCells-perPatient.ipynb).
- `adata.obs[[f"cNMF_{i}_score" for i in range(1,6)]+["sample_id","highlevel_refined"]].to_csv("/add/path/here/adata_cNMF_scores.csv")`: path to where to save the cNMF scores associated with each cell.
- `score_annotations.to_csv("/add/path/here/adata_cNMF_scores_wtop.csv")`: path to where to save the cNMF scores and the cNMF representative cells.
- `atac = sc.read_h5ad("/add/path/here/combined_atac.h5ad")`: path to .h5ad file containing the ATAC after Seurat analysis was run (see conversion from RDS to H5AD [here](https://github.com/vanallenlab/EAC-multiome/blob/main/code/R/scripts/preprocessing-snATAC/README.md)
- `adata_cnmf_scores = pd.read_csv("/add/path/here/adata_cNMF_scores.csv",index_col=0)`: path to the adata cNMF scores saved previously in the file.
- `most_corr_dir = pl.Path("/add/path/here/")`: path to where the most correlated ATAC regions associated with cNMF scores will be saved

##  8. `predict-score-from-atac.ipynb`

- `peak_info = pd.read_csv("/add/path/here/peaks_closestfeatures.csv").set_index("query_region")`:
- `atac = sc.read_h5ad("/add/path/here/combined_atac.h5ad")`: path to .h5ad file containing the ATAC after Seurat analysis was run (see conversion from RDS to H5AD [here](https://github.com/vanallenlab/EAC-multiome/blob/main/code/R/scripts/preprocessing-snATAC/README.md)
- scores = pd.read_csv("/add/path/here/adata_cNMF_scores_wtop.csv",index_col=0): path to the cNMF scores and representative cNMF cells saved in [`cohort-viz.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cohort-viz.ipynb).
- `most_corr_dir = pl.Path("/add/path/here/")`: path to where the most correlated ATAC regions associated with cNMF scores were saved in [`cohort-viz.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cohort-viz.ipynb).

#### RUN THE SCENIC PLUS ANALYSIS HERE
The next script analyses the results from the SCENIC+ analysis. More info on running the scripts [here](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/scripts/scenicplus/README.md).

## 9. `SCENICplus_analyze_cNMF.ipynb`

- `work_dir = pl.Path("/add/path/here")`:
- `refined_wcancer = pd.read_csv("/add/path/here/refined_annotations_wsampleid.csv",index_col=0)`: path to where the annotations were saved in [`cohort-viz.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cohort-viz.ipynb).
- `all_tfs = pd.read_csv("/add/path/here/DatabaseExtract_v_1.01.csv",index_col=0)`: path to where the Human Transcription Factor Database list of TFs was saved. (**DWNL**)
- `scores = pd.read_csv("/add/path/here/adata_cNMF_scores.csv",index_col=0)`: path to where the cNMF scores were saved in [`cohort-viz.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/cohort-viz.ipynb).
- `df.to_csv(f"/add/path/here/{state}_triad_corr.csv")`: path to where to save the results of the correlation analysis.
- `pd.DataFrame(df_toptfs).to_csv("/add/path/here/toptfs_top20.csv")`: path to where to save the top 20 candidate mTFs per program


