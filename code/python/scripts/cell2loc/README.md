# EAC-multiome

This subfolder contains the scripts to run the Cell2Location analysis as well as the bash scripts used to run them. The bash script are run on a SLURM system and are given as examples.

Make sure the environment you are in contains Cell2Location. You can create an environment from the provided [`cell2loc_env.yml`](https://github.com/vanallenlab/EAC-multiome/blob/main/cell2loc_env.yml) or create an environment with Cell2Location on your own. Data that needs to be downloaded prior to running the scripts is indicated by the (**DWNL**) tag.

**Note:** the [`Cell2Location.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/spatial-transcriptomics/1.%20Cell2Location.ipynb) notebook should be run before this step.

## `script-cell2loc.py`

- `spatial_dir = pl.Path("/add/path/here/")`: path to where the SpaceRanger output spatial data is saved. (**DWNL**). This should contain a folder for each sample with name "sample_name".
- `inf_aver = pd.read_csv("/add/path/here/mean-expression-per-celltype.csv",index_col=0,)`: path to the file with the mean expression per cell type estimated with the NB model in [`Cell2Location.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/spatial-transcriptomics/1.%20Cell2Location.ipynb). 
- `results_folder = f"/add/path/here/Cell2Location_results/{patient_name}"`:  path to where the Cell2Location results will be saved.