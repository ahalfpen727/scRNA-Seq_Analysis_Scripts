## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
grun.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE81nnn/GSE81076/suppl/GSE81076%5FD2%5F3%5F7%5F10%5F17%2Etxt%2Egz"))

## --------------------------------------------------------------------------
gse81076.df <- read.table(grun.fname, sep='\t', 
    header=TRUE, stringsAsFactors=FALSE, row.names=1)
dim(gse81076.df)

## --------------------------------------------------------------------------
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse81076.df))
table(donor.names)
plate.id <- sub("^D[0-9]+(.*)_.*", "\\1", colnames(gse81076.df))
table(plate.id)

## --------------------------------------------------------------------------
gene.symb <- gsub("__chr.*$", "", rownames(gse81076.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.symb[is.spike]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse81076.df <- gse81076.df[keep,]
rownames(gse81076.df) <- gene.ids[keep]
summary(keep)

## --------------------------------------------------------------------------
library(SingleCellExperiment)
sce.gse81076 <- SingleCellExperiment(list(counts=as.matrix(gse81076.df)),
	colData=DataFrame(Donor=donor.names, Plate=plate.id),
	rowData=DataFrame(Symbol=gene.symb[keep]))
isSpike(sce.gse81076, "ERCC") <- grepl("^ERCC-", rownames(gse81076.df)) 
sce.gse81076  

## --------------------------------------------------------------------------
library(scater)
sce.gse81076 <- calculateQCMetrics(sce.gse81076, compact=TRUE)
QC <- sce.gse81076$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
	HighSpike=sum(high.spike, na.rm=TRUE))

## --------------------------------------------------------------------------
discard <- low.lib | low.genes | high.spike
sce.gse81076 <- sce.gse81076[,!discard]
summary(discard)

## --------------------------------------------------------------------------
library(scran)
set.seed(1000)    
clusters <- quickCluster(sce.gse81076, method="igraph", min.mean=0.1)
table(clusters)
sce.gse81076 <- computeSumFactors(sce.gse81076, min.mean=0.1, clusters=clusters)
summary(sizeFactors(sce.gse81076))

## --------------------------------------------------------------------------
sce.gse81076 <- computeSpikeFactors(sce.gse81076, general.use=FALSE)
summary(sizeFactors(sce.gse81076, "ERCC"))

## --------------------------------------------------------------------------
sce.gse81076 <- normalize(sce.gse81076)

## ----var-gse81076, fig.cap="Variance of normalized log-expression values for each gene in the GSE81076 dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
block <- paste0(sce.gse81076$Plate, "_", sce.gse81076$Donor)
fit <- trendVar(sce.gse81076, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse81076, fit)

plot(dec$mean, dec$total, xlab="Mean log-expression", 
	ylab="Variance of log-expression", pch=16)
is.spike <- isSpike(sce.gse81076)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

## --------------------------------------------------------------------------
dec.gse81076 <- dec
dec.gse81076$Symbol <- rowData(sce.gse81076)$Symbol
dec.gse81076 <- dec.gse81076[order(dec.gse81076$bio, decreasing=TRUE),]
head(dec.gse81076)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(gse81076.df)
gc()

## --------------------------------------------------------------------------
muraro.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE85nnn/GSE85241/suppl",
    "GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz"))

## --------------------------------------------------------------------------
gse85241.df <- read.table(muraro.fname, sep='\t', 
    header=TRUE, row.names=1, stringsAsFactors=FALSE)
dim(gse85241.df)

## --------------------------------------------------------------------------
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse85241.df))
table(donor.names)
plate.id <- sub("^D[0-9]+\\.([0-9]+)_.*", "\\1", colnames(gse85241.df))
table(plate.id)

## --------------------------------------------------------------------------
gene.symb <- gsub("__chr.*$", "", rownames(gse85241.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.symb[is.spike]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse85241.df <- gse85241.df[keep,]
rownames(gse85241.df) <- gene.ids[keep]
summary(keep)

## --------------------------------------------------------------------------
sce.gse85241 <- SingleCellExperiment(list(counts=as.matrix(gse85241.df)),
	colData=DataFrame(Donor=donor.names, Plate=plate.id),
	rowData=DataFrame(Symbol=gene.symb[keep]))
isSpike(sce.gse85241, "ERCC") <- grepl("^ERCC-", rownames(gse85241.df)) 
sce.gse85241  

## --------------------------------------------------------------------------
sce.gse85241 <- calculateQCMetrics(sce.gse85241, compact=TRUE)
QC <- sce.gse85241$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
	HighSpike=sum(high.spike, na.rm=TRUE))

## --------------------------------------------------------------------------
discard <- low.lib | low.genes | high.spike
sce.gse85241 <- sce.gse85241[,!discard]
summary(discard)

## --------------------------------------------------------------------------
set.seed(1000)
clusters <- quickCluster(sce.gse85241, min.mean=0.1, method="igraph")
table(clusters)
sce.gse85241 <- computeSumFactors(sce.gse85241, min.mean=0.1, clusters=clusters)
summary(sizeFactors(sce.gse85241))
sce.gse85241 <- computeSpikeFactors(sce.gse85241, general.use=FALSE)
summary(sizeFactors(sce.gse85241, "ERCC"))
sce.gse85241 <- normalize(sce.gse85241)

## ----var-gse85241, fig.cap="Variance of normalized log-expression values for each gene in the GSE85241 dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
block <- paste0(sce.gse85241$Plate, "_", sce.gse85241$Donor)
fit <- trendVar(sce.gse85241, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse85241, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
	ylab="Variance of log-expression", pch=16)
is.spike <- isSpike(sce.gse85241)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

## --------------------------------------------------------------------------
dec.gse85241 <- dec
dec.gse85241$Symbol <- rowData(sce.gse85241)$Symbol
dec.gse85241 <- dec.gse85241[order(dec.gse85241$bio, decreasing=TRUE),]
head(dec.gse85241)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(gse85241.df)
gc()

## --------------------------------------------------------------------------
universe <- intersect(rownames(dec.gse85241), rownames(dec.gse81076))
mean.bio <- (dec.gse85241[universe,"bio"] + dec.gse81076[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

## --------------------------------------------------------------------------
rescaled <- multiBatchNorm(sce.gse85241[universe,], sce.gse81076[universe,])
rescaled.gse85241 <- rescaled[[1]]
rescaled.gse81076 <- rescaled[[2]]

## --------------------------------------------------------------------------
set.seed(100) 
original <- list(
    GSE81076=logcounts(rescaled.gse81076)[chosen,],
    GSE85241=logcounts(rescaled.gse85241)[chosen,]
)

# Slightly convoluted call to avoid re-writing code later.
# Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, approximate=TRUE)))
dim(mnn.out$corrected)

## --------------------------------------------------------------------------
mnn.out$batch

## --------------------------------------------------------------------------
mnn.out$pairs

## --------------------------------------------------------------------------
omat <- do.call(cbind, original)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- mnn.out$corrected
sce$Batch <- as.character(mnn.out$batch)
sce

## ----tsne-batch, fig.width=10, fig.asp=0.6, fig.cap="t-SNE plots of the pancreas datasets, before and after MNN correction. Each point represents a cell and is coloured by the batch of origin."----
set.seed(100)
# Using irlba to set up the t-SNE, for speed.
osce <- runPCA(sce, ntop=Inf, method="irlba")
osce <- runTSNE(osce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")

set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN")
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")

multiplot(ot, ct, cols=2)

## ----tsne-markers, fig.width=10, fig.height=10, fig.cap="t-SNE plots after MNN correction, where each point represents a cell and is coloured by its corrected expression of key marker genes for known cell types in the pancreas."----
ct.gcg <- plotTSNE(csce, colour_by="ENSG00000115263") + 
    ggtitle("Alpha cells (GCG)")
ct.ins <- plotTSNE(csce, colour_by="ENSG00000254647") + 
    ggtitle("Beta cells (INS)")
ct.sst <- plotTSNE(csce, colour_by="ENSG00000157005") + 
    ggtitle("Delta cells (SST)")
ct.ppy <- plotTSNE(csce, colour_by="ENSG00000108849") + 
    ggtitle("PP cells (PPY)")
multiplot(ct.gcg, ct.ins, ct.sst, ct.ppy, cols=2)

## --------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership, sce$Batch)

## ----tsne-cluster, fig.cap="t-SNE plot after MMN correction, where each point represents a cell and is coloured by its cluster identity."----
csce$Cluster <- factor(clusters$membership)
plotTSNE(csce, colour_by="Cluster")

## --------------------------------------------------------------------------
m.out <- findMarkers(sce, clusters$membership, block=sce$Batch,
    direction="up")        
demo <- m.out[["4"]] # looking at cluster 4 (probably alpha cells).
demo <- demo[demo$Top <= 5,]

library(org.Hs.eg.db)
data.frame(row.names=rownames(demo),
    Symbol=mapIds(org.Hs.eg.db, keytype="ENSEMBL", 
        keys=rownames(demo), column="SYMBOL"),
    Top=demo$Top, FDR=demo$FDR)

## ---- echo=FALSE, results="hide"-------------------------------------------
# Checking that cluster 4 is what we think it is.
stopifnot(all(sapply(demo["ENSG00000115263",-(1:3)], sign)==1))

## --------------------------------------------------------------------------
# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
design <- model.matrix(~sce$Batch)
design <- design[,-1,drop=FALSE]

m.alt <- findMarkers(sce, clusters$membership, design=design,
    direction="up")
demo <- m.alt[["4"]]
demo <- demo[demo$Top <= 5,]

data.frame(row.names=rownames(demo),
    Symbol=mapIds(org.Hs.eg.db, keytype="ENSEMBL", 
        keys=rownames(demo), column="SYMBOL"),
    Top=demo$Top, FDR=demo$FDR)

## --------------------------------------------------------------------------
sessionInfo()

