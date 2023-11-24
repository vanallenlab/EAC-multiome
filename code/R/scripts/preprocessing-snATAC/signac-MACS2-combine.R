library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(harmony)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)
library(TFBSTools)
library(patchwork)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

ct.annotations <- read.csv("add/path/here")

basedir <- "/add/path/here/"
macsdir <- "/add/path/here/"

sample.list <- c('Aguirre_EGSFR0148', 'Aguirre_EGSFR1732',
                 'Aguirre_EGSFR1938', 'Aguirre_EGSFR1982',
                 'Aguirre_EGSFR2218',
                 'Aguirre_EGSFR0074', 'CCG1153_4411', 
                 'CCG1153_4496262', 'CCG1153_6640539',
                 'Aguirre_EGSFR0128')

peaks.1 <- read.table(
  file = file.path(macsdir, sample.list[[1]] ,paste0(sample.list[[1]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.1 <- makeGRangesFromDataFrame(peaks.1)

peaks.2 <- read.table(
  file = file.path(macsdir, sample.list[[2]] ,paste0(sample.list[[2]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.2 <- makeGRangesFromDataFrame(peaks.2)

peaks.3 <- read.table(
  file = file.path(macsdir, sample.list[[3]] ,paste0(sample.list[[3]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.3 <- makeGRangesFromDataFrame(peaks.3)

peaks.4 <- read.table(
  file = file.path(macsdir, sample.list[[4]] ,paste0(sample.list[[4]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.4 <- makeGRangesFromDataFrame(peaks.4)

peaks.5 <- read.table(
  file = file.path(macsdir, sample.list[[5]] ,paste0(sample.list[[5]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.5 <- makeGRangesFromDataFrame(peaks.5)

peaks.6 <- read.table(
  file = file.path(macsdir, sample.list[[6]] ,paste0(sample.list[[6]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.6 <- makeGRangesFromDataFrame(peaks.3)

peaks.6 <- read.table(
  file = file.path(macsdir, sample.list[[6]] ,paste0(sample.list[[6]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.6 <- makeGRangesFromDataFrame(peaks.6)

peaks.7 <- read.table(
  file = file.path(macsdir, sample.list[[7]] ,paste0(sample.list[[7]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.7 <- makeGRangesFromDataFrame(peaks.7)

peaks.8 <- read.table(
  file = file.path(macsdir, sample.list[[8]] ,paste0(sample.list[[8]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.8 <- makeGRangesFromDataFrame(peaks.8)

peaks.9 <- read.table(
  file = file.path(macsdir, sample.list[[9]] ,paste0(sample.list[[9]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.9 <- makeGRangesFromDataFrame(peaks.9)

peaks.10 <- read.table(
  file = file.path(macsdir, sample.list[[10]] ,paste0(sample.list[[10]],".bed")),
  col.names = c("chr", "start", "end","name","score","xtra")
)
gr.10 <- makeGRangesFromDataFrame(peaks.10)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.1, gr.2, gr.3, gr.4, gr.5, gr.6, gr.7, gr.8, gr.9, gr.10))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

all.mds <- list()
for(sample.name in sample.list){
  all.mds[[sample.name]] <- read.table(
    file = file.path(basedir, sample.name, "per_barcode_metrics.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  all.mds[[sample.name]] <- all.mds[[sample.name]][all.mds[[sample.name]]$is_cell==1,]
  
}

all.samples <- list()
for(sample.name in sample.list){  
  frags <- CreateFragmentObject(
    path = file.path(basedir, sample.name, "atac_fragments.tsv.gz"),
    cells = rownames(all.mds[[sample.name]])
  )
  
  sample.counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    cells = rownames(all.mds[[sample.name]])
  )
  
  sample_assay <- CreateChromatinAssay(sample.counts, fragments = frags)
  sample <- CreateSeuratObject(sample_assay, assay = "ATAC", meta.data=all.mds[[sample.name]])
  
  subannot <- ct.annotations[which(ct.annotations$sample_id==sample.name),]
  subannot$X <- substr(subannot$X, 1, nchar(subannot$X)-2)
  rownames(subannot) <- subannot$X
  sample <- AddMetaData(object=sample, metadata=subannot)
  sample$dataset <- sample.name
  
  all.samples[[sample.name]] <- sample
}

respath <- "/add/path/here/"

for(sample.name in sample.list){
  print(sample.name)
  saveRDS(object = all.samples[[sample.name]], file=file.path(respath, "sample-level-atac-rds", paste0(sample.name, ".rds")))
}

all.samples <- list()
for(sample.name in sample.list){
  print(sample.name)
  all.samples[[sample.name]] <- readRDS(file=file.path(respath, "sample-level-atac-rds", paste0(sample.name, ".rds")))
}

combined <- merge(
  x = all.samples[[sample.list[[1]]]],
  y = c(all.samples[[sample.list[[2]]]],
        all.samples[[sample.list[[3]]]],
        all.samples[[sample.list[[4]]]],
        all.samples[[sample.list[[5]]]],
        all.samples[[sample.list[[6]]]],
        all.samples[[sample.list[[7]]]],
        all.samples[[sample.list[[8]]]],
        all.samples[[sample.list[[9]]]],
        all.samples[[sample.list[[10]]]]),
  add.cell.ids = c("1","2","3","4","5","6","7","8","9","10")
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

Annotation(combined) <- annotations

combined <- NucleosomeSignal(combined)
combined <- TSSEnrichment(combined)

VlnPlot(
  object = combined,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


combined <- subset(
  x = combined,
  subset = nCount_ATAC < 500000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
p2 <- DimPlot(combined, group.by = 'highlevel_annotation', pt.size = 0.1)
p1 | p2

combined <- RegionStats(object=combined, genome = BSgenome.Hsapiens.UCSC.hg38, assay="ATAC")
saveRDS(object = combined, file = file.path(respath, "combined_atac.rds"))

subset_cells <- colnames(combined)[which(combined@meta.data$refined_wcancer %in% c("cNMF_1","cNMF_3","cNMF_4","Carcinoma_undefined"))]
malcombined <- subset(x = combined, cells=subset_cells)

malcombined <- RunTFIDF(malcombined)
malcombined <- FindTopFeatures(malcombined, min.cutoff = 20)
malcombined <- RunSVD(malcombined)
malcombined <- RunUMAP(malcombined, dims = 2:50, reduction = 'lsi')

saveRDS(object = malcombined, file = file.path(respath, "malignant_combined_atac.rds"))

p1 <- DimPlot(malcombined, group.by = 'dataset', pt.size = 0.1)
p2 <- DimPlot(malcombined, group.by = 'refined_wcancer', pt.size = 0.1)
p1 | p2

malcombined <- readRDS(file = file.path(respath, "malignant_combined_atac.rds"))

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
malcombined <- AddMotifs(
  object = malcombined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#### WILL BE USEFUL FOR SEACELL COMPUTATION
# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(malcombined), sep = c(":", "-")),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  score = TRUE,
  sep = c(":", "-"),
  use.counts = FALSE
)
motif_pwm = as.matrix(motif.matrix)
write.csv(motif_pwm, "/add/path/here/motif_pwm_score.csv")


Idents(malcombined) <- "refined_wcancer"

diffpeakdir <- "/add/path/here"

diffcnmfs <- c("cNMF_1","cNMF_3","cNMF_4")

for(id.1 in diffcnmfs){
  da_peaks <- FindMarkers(
    object = malcombined,
    ident.1 = id.1,
    ident.2 = setdiff(diffcnmfs, c(id.1)),
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
  )
  
  write.csv(da_peaks, file=file.path(diffpeakdir,paste0(id.1,".csv")))
  
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
  
  enriched.motifs <- FindMotifs(
    object = malcombined,
    features = top.da.peak
  )
  
  write.csv(enriched.motifs, file=file.path(diffpeakdir,paste0(id.1,"_enrichedmotifs.csv")))
}

for(id.1 in diffcnmfs){
  da_peaks <- read.csv(file.path(diffpeakdir,paste0(id.1,".csv")), row.names = 1)
  openpeaks <- rownames(da_peaks[da_peaks$p_val_adj <0.05 & da_peaks$avg_log2FC > 1.5, ])
  closestfeat <- ClosestFeature(malcombined, openpeaks)
  write.csv(closestfeat, file=file.path(diffpeakdir,paste0(id.1,"_closestfeatures.csv")))
}

openpeaks <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
closestfeat <- ClosestFeature(malcombined, openpeaks)
CoveragePlot(
  object = malcombined,
  region = c("RAP1B", "IRAK3","PTPRB"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  ncol = 1
)

peaks = CallPeaks(
  object = malcombined,
  macs2.path = "/Users/josephineyates/opt/anaconda3/envs/spatial/bin/macs2",
  group.by = 'refined_wcancer')

saveRDS(object = peaks, file = glue('data/single_sample_processing/{sample_id}/objects/{sample_id}_clusterlevel_peaks.rds'))