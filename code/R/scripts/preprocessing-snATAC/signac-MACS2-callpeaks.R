library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)


plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

basedir <- "/add/path/here"
macsdir <- "/add/path/here"

sample.list <- c('Aguirre_EGSFR0148', 'Aguirre_EGSFR1732',
                 'Aguirre_EGSFR1938', 'Aguirre_EGSFR1982',
                 'Aguirre_EGSFR2218',
                 'Aguirre_EGSFR0074', 'CCG1153_4411', 
                 'CCG1153_4496262', 'CCG1153_6640539',
                 'Aguirre_EGSFR0128')


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

for(sample.name in sample.list){
  print(sample.name)
  dir.create(file.path(macsdir, sample.name), showWarnings = FALSE)
  
  md <- read.table(
    file = file.path(basedir, sample.name, "per_barcode_metrics.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  md <- md[md$is_cell==1,]
  
  counts <- Read10X_h5(paste0("/add/path/here",sample.name,"/filtered_feature_bc_matrix.h5"))
  fragpath <- paste0("/add/path/here",sample.name,"/atac_fragments.tsv.gz")

  
  # create ATAC assay and add it to the object
  sample_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  sample.obj <- CreateSeuratObject(sample_assay, assay = "ATAC", meta.data=md)
  
  # call peaks using MACS2
  peaks <- CallPeaks(sample.obj, 
                     macs2.path="/add/path/here/macs2", )
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  
  export(peaks, file.path(macsdir, sample.name,paste0(sample.name,".bed")), "bed")

}