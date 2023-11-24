# EAC-multiome

## `signac-MACS2-callpeaks.R`
- `basedir <- "/add/path/here"`
- `macsdir <- "/add/path/here"`
- `counts <- Read10X_h5(paste0("/add/path/here",sample.name,"/filtered_feature_bc_matrix.h5"))`
- `fragpath <- paste0("/add/path/here",sample.name,"/atac_fragments.tsv.gz")`
- `macs2.path="/add/path/here/macs2"`

## `signac-MACS2-combine.R`
- `ct.annotations <- read.csv("add/path/here")`
- `basedir <- "/add/path/here/"`
- `macsdir <- "/add/path/here/"`
- `respath <- "/add/path/here/"`
- `write.csv(motif_pwm, "/add/path/here/motif_pwm_score.csv")`
- `diffpeakdir <- "/add/path/here"`
