# EAC-multiome

This subfolder contains the scripts to run the pyscenicplus analysis as well as the bash scripts used to run them. The bash script are run on a SLURM system and are given as examples.

First, SCENIC+ must be installed using the instructions [here](https://github.com/aertslab/scenicplus). 
Then, the auxiliary resources needed for running SCENIC+ must be downloaded. All the information about files that need to be downloaded for SCENIC and SCENIC+ is given in the [main README](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md). Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

**Note:** the [`signac-MACS2-callpeaks.R`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/R/scripts/preprocessing-snATAC/signac-MACS2-callpeaks.R) should be run before this step, as well as all python notebooks from the preprocessing folder and the analysis folder up until the 8. step.

## 1. `run-pre-pyscenic-script.py`

- `work_dir = pl.Path("/add/path/here")`: the path to where the results will be stored and then used for the SCENIC+ analysis script [`9. SCENICplus-analyze-cNMF.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/9.%20SCENICplus-analyze-cNMF.ipynb)
- `atacdir = pl.Path("/add/path/here/")`: path to where the output from the CellRanger ARC pipeline is kept, in particular where files linked to ATAC are. This should contain a folder for each sample that contains at least `atac_fragments.tsv.gz` and `per_barcode_metrics.csv`. (**DWNL**)
- `macsdir = pl.Path("/add/path/here/")`: path to where the MACS2 output from MACS2 peak calling is stored after the script [`signac-MACS2-callpeaks.R`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/R/scripts/preprocessing-snATAC/signac-MACS2-callpeaks.R). This should contain a folder for each sample containing the `.bed` files for the MACS2 called peaks. 
- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/2.%20tme-cleaning-analysis.ipynb)
- `path_to_blacklist = "/add/path/here/hg38-blacklist.v2.bed"`: path to where the hg38 blacklisted regions `.bed` is located. (**DWNL**)
- `tmp_dir="/add/path/here`: path to where the temporary data (e.g., RAY spills) can be saved. This needs to have a decent amount of memory, beware!
- `auxiliarypath = pl.Path("/add/path/here/")`: path to the folder where the SCENIC+ resources are saved (**DWNL**)

## `run-scenicplus.py`

- `work_dir = pl.Path("/add/path/here/")`: path to where the results have been saved up until now in the script [`1. run-pre-pyscenic-script.py`](TBD)
- `tmp_dir="/add/path/here`: path to where the temporary data (e.g., RAY spills) can be saved. This needs to have a decent amount of memory, beware!
- `annot_dir = pl.Path("/add/path/here/")`: path to the folder where the annotations are saved to be able to run the searchspace. (**DWNL**)
- `adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")`: path to where the full cleaned cohort `.h5ad` object was saved in [`tme-cleaning-analysis.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/analysis/2.%20tme-cleaning-analysis.ipynb)
- `tf_file = "/add/path/here/utoronto_human_tfs_v_1.01.txt"`: path to where the file with the full list of TFs created by Lambert et al. is located (**DWNL**)