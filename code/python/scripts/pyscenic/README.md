# EAC-multiome

This subfolder contains the instructions to run the SCENIC analysis on external datasets. All instructions to download necessary files are located [here](https://github.com/vanallenlab/EAC-multiome/blob/main/README.md)

SCENIC was run using the command line tool on a HPC. Here we provide the CLI instructions to get the results.

The study name is indicated by $STUDY. Replace it by Carroll, Croft, or Luo depending on the run.

Path to the files need to be replaced where `/add/path/here` is indicated.
First run 
```
arboreto_with_multiprocessing.py \
    /add/path/here/adata_$STUDY_filtered_scenic.loom \
    /add/path/here/utoronto_human_tfs_v_1.01.txt \
    --method grnboost2 \
    --output adj_$STUDY.tsv \
    --num_workers 32 \
    --seed 777
```

`/add/path/here/adata_$STUDY_filtered_scenic.loom`corresponds to the loom file saved in the validation script [`5. runpySCENIC-external part 1.ipynb`](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/validation/5.%20runpySCENIC-external%20part%201.ipynb) for the study in question (Carroll, Luo, or Croft) and `/add/path/here/utoronto_human_tfs_v_1.01.txt` corresponds to the list of human TFs curated by Lambert et al. 2018, file 'List of Lambert et al. TF names'.

Then run 
```
pyscenic ctx adj_$STUDY.tsv \
    /add/path/here/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /add/path/here/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /add/path/here/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname adata_$STUDY_filtered_scenic.loom \
    --output reg_$STUDY.csv \
    --mask_dropouts \
    --num_workers 32
```
` /add/path/here/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` corresponds to the ranking of motifs, file 'Screen v10 hg38 database, SCENIC, #1'; `/add/path/here/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`corresponds to the ranking of motifs, file 'Screen v10 hg38 database, SCENIC, #2'; `/add/path/here/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`corresponds to the motif annotations, file 'Motif v10 annotation, SCENIC+'. 

Finally run 
```
pyscenic aucell \
    adata_$STUDY_filtered_scenic.loom \
    reg_$STUDY.csv \
    --output pyscenic_$STUDY_output.loom\
    --num_workers 32
```
