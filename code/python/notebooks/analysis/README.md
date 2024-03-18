# EAC-multiome

This folder contains all notebooks used for the analysis of the discovery cohort. 

The information to download the RNA counts for the datasets as well as all auxiliary files are located in the [related subsections of the README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

These notebooks should be run in the following order:
## `cohort_plot.ipynb`

- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**)
- `toy_comut.figure.savefig("add/path/here/full_cohort_plot.svg", dpi=200)`: path to where to save the Figure.

## `tme-cleaning-analysis.ipynb`

- `datadir = pl.Path("add/path/here/celltype_h5ad_files/")`: path to the folder that contains the output of the preprocessing notebooks. Should contain one `.h5ad` file per patient.
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**)
- `refined_annotations.to_csv("add/path/here/refined_annotations.csv")`: path to where to save the refined annotations derived from this notebook. 
- `adata.write("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object will be saved. 

## `plotting_TME_subsets.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/tme-cleaning-analysis.ipynb)
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`: path to the clinical data of the discovery cohort (**DWNL**).


## `cNMFCancerCells-perPatient.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`
- `program_dir = pl.Path("/add/path/here/cNMF_malignant_per_patient/")`
- `marker_genes[cl].to_csv(f"/add/path/here/cNMF_{cl}.csv")`
- `itay_MPs = pd.read_csv("/add/path/here/ItayTiroshHeterogeneityMPs.csv")`
- `gene_sets='/add/path/here/h.all.v7.4.symbols.gmt',`

## `cohort_viz.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv",index_col=0)`
- `refined_wcancer = pd.read_csv("/add/path/here/refined_wCNMF_programs_and_sampleid.csv",index_col=0)`
- `marker_genes[cl] = pd.read_csv(f"/add/path/here/cNMF_{cl}.csv",index_col=0)`
- `atac = sc.read_h5ad("/add/path/here/combined_atac.h5ad")`
- `DAR_res_dir = pl.Path("/add/path/here")`
- `cnmf_res_dir = pl.Path("/add/path/here")`

## `investigate_mTF_states.ipynb`
- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv",index_col=0)`
- `refined_wcancer = pd.read_csv("/add/path/here/refined_wCNMF_programs_and_sampleid.csv",index_col=0)`
- `marker_genes[cl] = pd.read_csv(f"/add/path/here/cNMF_{cl}.csv",index_col=0)`

## `run-milopy.ipynb`

- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`
- `clinical = pd.read_csv("/add/path/here/EAC_clinical_info.csv",index_col=0)`
- `refined_annotations = pd.read_csv("/add/path/here/refined_annotations.csv",index_col=0)`
- `refined_wcancer = pd.read_csv("/add/path/here/refined_wCNMF_programs_and_sampleid.csv",index_col=0)`

## `SCENICplus_analyze_cNMF.ipynb`

- `work_dir = pl.Path("/add/path/here")`

## `score_celllines.ipynb`

- `celline_dir = pl.Path("/add/path/here/")`
- `signature_dir = pl.Path("/add/path/here")`

## `ScoreEACinTCGA.ipynb`

- `tcga_dir = pl.Path("/add/path/here")`
- `eac_tcga_dir = pl.Path("/add/path/here")`
- `signature_dir = pl.Path("/add/path/here")`
- `path_reports = pd.read_csv("/add/path/here/TCGA_Reports.csv",index_col=0)`
- `pni_info.to_csv("/add/path/here/PNI_status_annotated_manually.csv")`
- `pni_info = pd.read_csv("/add/path/here/PNI_status_annotated_manually.csv",index_col=0)`
- `bp_deconvolved[cancer] = pd.read_csv(f"/add/path/here/{name_mapping[cancer]}_tumor_component_counts.csv",index_col=0)`
- `bp_purity[cancer] = pd.read_csv(f"/add/path/here/{name_mapping[cancer]}_purity.csv",index_col=0)`
- `gencode = pd.read_csv("/add/path/here/gencode_annot_length.csv",index_col=0).set_index("gene_name")`

## `TCGAlinkCNAburden.ipynb`

- `tcga_dir = pl.Path("/add/path/here")`
- `signature_dir = pl.Path("/add/path/here")`


