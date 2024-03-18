library(BayesPrism)
library(Seurat)

path <- "/add/path/here"

expr <- "/add/path/here/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt"

seuratObj <- readRDS(paste0(path, "add/path/here/sce.pub.Rds"))

data <- read.table(paste0(path, expr), sep='\t', row.names=1, header=TRUE) 

#prepare data for bayesprism:

bk.dat <- t(data)
bk.dat <- bk.dat[,2:18845]

### EAC
celltypes = seuratObj$subcompartment
cellstates = seuratObj$cellsubtype
idx_escc <- grep("ESCC-", cellstates)
idx_eac <- grep("EAC-", cellstates)
cellstates[c(idx_eac)] <- "EAC"
cellstates[c(idx_escc)] <- "ESCC"
celltypes[c(idx_eac,idx_escc)] <- "tumor"

### EAC
sc.dat <- as.matrix(t(seuratObj@assays$RNA@counts))
bk.dat <- bk.dat[,!duplicated(colnames(bk.dat))] #remove duplicated genes

#create BP object:

myPrism <- new.prism(reference=sc.dat, mixture=bk.dat, input.type="count.matrix",
                     cell.type.labels = celltypes, cell.state.labels = cellstates,
                     key="tumor", outlier.cut = 0.001, outlier.fraction=0.1)

#run BP:

bp.res <- run.prism(prism=myPrism, n.cores=15, update.gibbs = TRUE)


#save the object for later:
save(bp.res, file=paste0(path, "/add/path/here/eac-gse207526.bp.res.RData"))

