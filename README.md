# EAC-multiome

This repository accompanies the work "[TBD]" (ref)
It contains all the code necessary to reproduce the analyses. 
Each subsection contains a README that contains a description of the path placeholders. 

## How to reproduce the analyses

To reproduce the analyses, one first needs to download the data available at XXX ([Instructions](link))

Then, one needs to run the scripts in order, as some intermediate files generated by the scripts will be re-used in subsequent scripts. 

For the snRNA-seq analysis, one needs to 
- run the preprocessing for each patient ([Instructions](link))
- run the preprocessing scripts for the external dataset to put them in the right format ([Instructions](link))
- run the analyses ([Instructions](link))
- run the pyscenic plus scripts ([Instructions](link))
- run the pyscenic scripts on external datasets ([Instructions](link))
- run the scripts for external validation ([Instructions](link))

For the snATAC-seq analysis, one needs to 
- run the prepocessing scripts ([Instructions](link))
- run BayesPrism for TCGA validation ([Instructions](link))

## Patient ID to sample ID mapping

In the original paper, for simplicity patients are referred to as P1 through P10. In the scripts/notebooks the patients are referred to using their sample ID. The mapping is provided below.
| Patient ID             | Sample ID                  |
|---------------------|--------------------------------|
| P1   | CCG1153_4496262   |
| P2   | CCG1153_6640539   |
| P3   | CCG1153_4411   |
| P4   | Aguirre_EGSFR0074   |
| P5   | Aguirre_EGSFR0148   |
| P6   | Aguirre_EGSFR1732   |
| P7   | Aguirre_EGSFR0128   |
| P8   | Aguirre_EGSFR1938   |
| P9   | Aguirre_EGSFR1982   |
| P10   | Aguirre_EGSFR2218   |

## Links to data used in the study

| Dataset             | Link to paper                  | Link to download                 | Remarks          | 
|---------------------|--------------------------------|---------------------------------|------------------|
| Discovery dataset    | [Yates et al., ???](TBD)    | [Download](http://example.com)  |     |
| Single-cell, Carroll et al.     | [Carroll et al., 2023](https://www.sciencedirect.com/science/article/pii/S1535610823002167?via%3Dihub)    | [Download](https://ega-archive.org/datasets/EGAD00001009401)  | Need to request access to data through EGA  / contact author (thomas.carroll@alumni.rice.edu)  |
| Bulk, Carroll et al., RNA   | [Carroll et al., 2023](https://www.sciencedirect.com/science/article/pii/S1535610823002167?via%3Dihub)    | [Download](https://ega-archive.org/datasets/EGAD00001009399)  | Need to request access to data through EGA  / contact author (thomas.carroll@alumni.rice.edu)  |
| Bulk, Carroll et al., Clinical   | [Carroll et al., 2023](https://www.sciencedirect.com/science/article/pii/S1535610823002167?via%3Dihub)     | [Download](https://bitbucket.org/licroxford/carroll_etal_2023/src/master/supplementary_files/Table_S8_papermetadata.xlsx)  | Inoperable cohort info is located [here](https://bitbucket.org/licroxford/carroll_etal_2023/src/master/supplementary_files/SuppInfo_OperablePts_LUD2015scRNAseq.xlsx)   |
| Single-cell, Croft et al.     | [Croft et al., 2022](https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-022-01666-x)    | [Download](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222078)  | Need to request single-cell annotations from author (w.d.croft@bham.ac.uk)    |
| Bulk, Hoefnagel et al., RNA     | [Hoefnagel et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9496882/)    | [Download](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207527)  |     |
| Bulk, Hoefnagel et al., Clinical     | [Hoefnagel et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9496882/)    | NA | Need to request from the author (sanne_hoefnagel@live.nl) |
| Bulk, TCGA, RNA (FPKM)| [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)    | [Download](https://xenabrowser.net/datapages/?dataset=TCGA-ESCA.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)  |  Used in the general TCGA analysis script, file named "TCGA-ESCA.htseq_fpkm-uq.tsv.gz" |
| Bulk, TCGA, RNA (Raw counts)| [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)    | [Download](https://xenabrowser.net/datapages/?dataset=TCGA-ESCA.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)  |  Used as a basis to deconvolve for BayesPrism  |
| Bulk, TCGA, Clinical #1     | [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)    | [Download](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.ESCA.sampleMap%2FESCA_clinicalMatrix)  | This is the general clinical+phenotypical info, named "TCGA.ESCA.sampleMap_ESCA_clinicalMatrix"  |
| Bulk, TCGA, Clinical #2     | [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)  | [Download](https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM91_ESM.xlsx?_gl=1*mw3dke*_ga*OTAwNDk3MDU2LjE3MDg2ODk0Mzk.*_ga_B3E4QL2TPR*MTcxMDc3MDI1My4yMS4wLjE3MTA3NzAzNTcuMC4wLjA.)  | This is the clinical info provided in the original paper, need to save as "ESCA_Nature_clinicalinfo.csv"   |
| Bulk, TCGA, Clinical #3     | [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)  | [Download](https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp)  | This is the survival information, file named "Survival_SupplementalTable_S1_20171025_xena_sp" |
| Bulk, TCGA, Clinical #4     | [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)  | [Download](https://api.gdc.cancer.gov/data/66dd07d7-6366-4774-83c3-5ad1e22b177e)  | This is the HRD information, file named "TCGA.HRD_withSampleID.txt" |
| Bulk, TCGA, ABSOLUTE purity  | [The Cancer Genome Atlas Research Network, 2017](https://www.nature.com/articles/nature20805)  | [Download](https://api.gdc.cancer.gov/data/66dd07d7-6366-4774-83c3-5ad1e22b177e)  | This is the ABSOLUE-estimate purity used for assessment of BayesPrism deconvolution, file named "TCGA_absolute_purity.txt" in the script |
| Single-cell, Luo et al.     | [Luo et al., 2022](https://www.nature.com/articles/s41467-022-34395-2)    | [Download](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210347)  | Need to download counts and metadata at the same time from this link    |


## Links to other data to download
| Data             | Needed for what script?                  | Description               | Link to download          | 
|---------------------|--------------------------------|---------------------------------|------------------|
| Gene Mapping    | R/scripts/BayesPrism/runBPrism.R    | Gene probe map fro the UCSC Xena browser that maps ENCODE to official gene ID |   [Download](https://github.com/ucscXena/xena-GDC-ETL/blob/master/xena_gdc_etl/resources/gencode.v22.annotation.gene.probeMap)  |
| GENCODE annotations   | python/notebooks/preprocessing-snRNA/XXXX.ipynb (where XXX is any sample name) | Subset of gencode annotations v41 |   [Download](TBD) or [link to original GTF file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz)  |
| Gene Programs from Gavish et al.   | python/notebooks/analysis/5. cNMFCancerCells-perPatient.ipynb | Signature genes derived in the [Gavish et al. paper](https://www.nature.com/articles/s41586-023-06130-4) |   [Download](TBD) or [link to original Excel file](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06130-4/MediaObjects/41586_2023_6130_MOESM6_ESM.xlsx), the .csv corresponds to the first sheet only |
| MSigDB Hallmarks of cancer GMT  | python/notebooks/analysis/5. cNMFCancerCells-perPatient.ipynb | This is the file to run GSEA on the hallmarks of cancer|   [Download](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/h.all.v7.4.symbols.gmt) |
| List of human transcription factors | python/notebooks/analysis/9. SCENICplus-analyze-cNMF.ipynb | This file contains all known human transcription factors as defined in the [Lambert et al. paper](https://www.sciencedirect.com/science/article/pii/S0092867418301065)  |   [Download](http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv) |
| Cell cycle genes | python/notebooks/validation/3. Carroll-validation-set.ipynb | This file contains cell cycle genes used by Scanpy |   [Download](https://github.com/scverse/scanpy_usage/blob/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt) |
| Marker genes of Barrett's esophagus cell types | python/notebooks/validation/4. compare-Nowicki-BE.ipynb | Set of fi!les, each containing marker genes of the Barrett's esophagus non-immune or stromal cell types |   [Download](TBD) or [Original paper tables](https://www.science.org/doi/suppl/10.1126/science.abd1449/suppl_file/science.abd1449_tables_s1_to_s12.zip); marker genes are derived from Suppl Table 7 |


