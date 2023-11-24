# EAC-multiome

## `CAFatlas-fibroblast-validation.ipynb`

- `metadata = pd.read_csv("/add/path/here/GSE210347_meta.txt.gz", sep="\t",index_col=0)`
- `adata = sc.read_h5ad("/add/path/here/GSE210347_counts.h5ad")`
- `fibadata.write_h5ad("/add/path/here/GSE210347_fibroblast_counts.h5ad")`

## `Carroll_validation_set.ipynb`

- `cell_cycle_genes = [x.strip() for x in open('/add/path/here/regev_lab_cell_cycle_genes.txt')]`
- `adata = sc.read_h5ad("/add/path/here/Carroll_singlecell/Carroll_EAC_raw.h5ad")`
- `signature_dir = pl.Path("/add/path/here")`
- `clinical = pd.read_csv("/add/path/here/carroll_clinical.csv", index_col=0)`
- `gencode = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`
- `bulk = pd.read_csv("/add/path/here/Carroll_singlecell/LUD2015-005_RNAseq_featureCounts.tsv",sep="\t",index_col=0)`

## `external_GSE207527.ipynb`

- `gex_df = pd.read_csv("/add/path/here/GSE207526/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt",sep="\t").iloc[1:,:].T`
- `gencode_mapping = pd.read_csv("/add/path/here/gencode_v41_positions.csv",index_col=0)`
- `gse = GEOparse.get_GEO(geo="GSE207526", destdir="/add/path/here")`
- `survival_clin = pd.read_csv("/add/path/here/data.SPSS.subselect.txt", sep="\t", index_col=0).set_index("FileName.GenomeScan")`
- `signature_dir = pl.Path("/add_path_here/")`

## `Croft_validation_set.ipynb`

- `datadir = pl.Path("/add/path/here")`
- `adata.write("/add/path/here/GSE222078_adata.h5ad")`
- `signature_dir = pl.Path("/add/path/here")`
- `adata = sc.read_h5ad("/add/path/here/GSE222078_adata.h5ad")`

