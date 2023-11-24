library(BayesPrism)
library(Seurat)
#library(SeuratDisk)

#source TCGA: https://xenabrowser.net/datapages/?dataset=TCGA-ESCA.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

geneMap <- "/add/path/here"
path <- "add/path/here"

expr <- "auxiliary_data/tcga/TCGA-ESCA.htseq_counts.tsv"

seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))

tcga <- read.table(paste0(path, expr), sep='\t', row.names=1, header=TRUE) #it is log2(count+1)

tcga <- (2^tcga) -1


#prepare data for bayesprism:

bk.dat = t(tcga[1:60483,]) # to omit info about mapped reads etc

### EAC
celltypes = seuratObj$subcompartment
cellstates = seuratObj$cellsubtype

idx_escc <- grep("ESCC-", cellstates)
idx_eac <- grep("EAC-", cellstates)
celltypes[c(idx_eac,idx_escc)] <- "tumor"
cellstates[c(idx_eac)] <- "EAC"
cellstates[c(idx_escc)] <- "ESCC"

sc.dat <- as.matrix(t(seuratObj@assays$RNA@counts))

gene2id <- read.table(paste0(path, geneMap), sep='\t', row.names=1, header=TRUE)

colnames(bk.dat) <- gene2id[colnames(bk.dat),"gene"] #to make sure bk and sc data have the same gene names

bk.dat <- bk.dat[,!duplicated(colnames(bk.dat))] #remove duplicated genes

#create BP object:

myPrism <- new.prism(reference=sc.dat, mixture=bk.dat, input.type="count.matrix",
                     cell.type.labels = celltypes, cell.state.labels = cellstates,
                     key="tumor", outlier.cut = 0.001, outlier.fraction=0.1)


#run BP:

bp.res <- run.prism(prism=myPrism, n.cores=15, update.gibbs = TRUE)


#save the object for later:
save(bp.res, file=paste0(path, "deconv_BayesPrism/eac.bp.res.RData"))

