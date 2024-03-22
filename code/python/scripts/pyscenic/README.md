# EAC-multiome

This subfolder contains the instructions to run the SCENIC analysis on external datasets.

SCENIC was run using the command line tool on a HPC. Here we provide the CLI instructions to get the results.

### 1. Croft et al.

Path to the files need to be replaced where `/add/path/here` is indicated.
First run 
```
arboreto_with_multiprocessing.py \
    /add/path/here/adata_croft_filtered_scenic.loom \
    /add/path/here/utoronto_human_tfs_v_1.01.txt \
    --method grnboost2 \
    --output adj_croft.tsv \
    --num_workers 32 \
    --seed 777
```

`/add/path/here/adata_croft_filtered_scenic.loom`corresponds to the loom file saved in the validation script 