## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
base.path <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2834nnn/GSM2834500/suppl"
barcode.fname <- bfcrpath(bfc, file.path(base.path, 
    "GSM2834500%5FG%5F1%5Fbarcodes%2Etsv%2Egz"))
gene.fname <- bfcrpath(bfc, file.path(base.path,
    "GSM2834500%5FG%5F1%5Fgenes%2Etsv%2Egz"))
counts.fname <- bfcrpath(bfc, file.path(base.path,
    "GSM2834500%5FG%5F1%5Fmatrix%2Emtx%2Egz"))

## --------------------------------------------------------------------------
library(scater)
library(Matrix)
gene.info <- read.table(gene.fname, stringsAsFactors=FALSE)
colnames(gene.info) <- c("Ensembl", "Symbol")
sce <- SingleCellExperiment(
    list(counts=as(readMM(counts.fname), "dgCMatrix")), 
    rowData=gene.info, 
    colData=DataFrame(Barcode=readLines(barcode.fname))
)

## --------------------------------------------------------------------------
rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$Ensembl, rowData(sce)$Symbol)
colnames(sce) <- sce$Barcode
sce

## --------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
chrloc <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keytype="GENEID", 
    keys=rowData(sce)$Ensembl, column="CDSCHROM")
rowData(sce)$Chr <- chrloc

## --------------------------------------------------------------------------
is.mito <- rowData(sce)$Chr == "chrM"
summary(is.mito)
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(is.mito)))

## --------------------------------------------------------------------------
low.lib <- isOutlier(sce$total_counts, log=TRUE, nmads=3, type="lower")
low.nexprs <- isOutlier(sce$total_features_by_counts, log=TRUE, nmads=3, type="lower")
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
discard <- low.lib | low.nexprs | high.mito
DataFrame(LowLib=sum(low.lib), LowNum=sum(low.nexprs), HighMito=sum(high.mito), 
    Discard=sum(discard), Kept=sum(!discard))

## --------------------------------------------------------------------------
sce <- sce[,!discard]

## --------------------------------------------------------------------------
library(scran)
set.seed(1000)
clusters <- quickCluster(sce, method="igraph", min.mean=0.1)
table(clusters)
sce <- computeSumFactors(sce, clusters=clusters, min.mean=0.1)
summary(sizeFactors(sce))

## --------------------------------------------------------------------------
sce <- normalize(sce)
assayNames(sce)

## ----varplot, fig.cap="Variance of the log-expression values as a function of the mean log-expression in the mammary gland data set. Each point represents a gene, and the red line corresponds to Poisson variance."----
tech.trend <- makeTechTrend(x=sce)
fit <- trendVar(sce, use.spikes=FALSE)
plot(fit$mean, fit$var, pch=16, 
    xlab="Mean log-expression",
    ylab="Variance of log-expression")
curve(tech.trend(x), add=TRUE, col="red")

## --------------------------------------------------------------------------
set.seed(12345)
sce <- denoisePCA(sce, technical=tech.trend, approximate=TRUE)
ncol(reducedDim(sce))

## --------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="PCA", k=25)
sce$Cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce$Cluster)

## ----tsneclust, fig.cap="t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its assigned cluster identity."----
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA")
plotTSNE(sce, colour_by="Cluster")

## --------------------------------------------------------------------------
dbl.out <- doubletCluster(sce, sce$Cluster)
dbl.out

## ----heatclust, fig.cap="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 7 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend."----
markers <- findMarkers(sce, sce$Cluster, direction="up")
dbl.markers <- markers[["7"]]
chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(sce, columns=order(sce$Cluster), colour_columns_by="Cluster", 
    features=chosen, cluster_cols=FALSE, center=TRUE, symmetric=TRUE, 
    zlim=c(-5, 5), show_colnames=FALSE)

## ---- echo=FALSE, results="hide"-------------------------------------------
# Checking that we've picked the correct cluster.
acta2 <- sapply(dbl.markers["Acta2", -(1:3)], sign)
csn2 <- sapply(dbl.markers["Csn2", -(1:3)], sign)

below <- acta2 < 0
stopifnot(all(csn2[below] == 1))
below <- csn2 < 0
stopifnot(all(acta2[below] == 1))

## ----markerexprs, fig.asp=0.5, fig.width=10, fig.cap="Distribution of log-normalized expression values for _Acta2_ and _Csn2_ in each cluster. Each point represents a cell."----
plotExpression(sce, features=c("Acta2", "Csn2"), 
    x="Cluster", colour_by="Cluster")

## --------------------------------------------------------------------------
set.seed(100)
dbl.dens <- doubletCells(sce, approximate=TRUE)
summary(dbl.dens)

## ----denstsne, fig.cap="t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its doublet density."----
sce$DoubletScore <- dbl.dens
plotTSNE(sce, colour_by="DoubletScore")

## ----densclust, fig.cap="Distribution of doublet scores for each cluster in the mammary gland data set. Each point is a cell."----
plotColData(sce, x="Cluster", y="DoubletScore", colour_by="Cluster")

## --------------------------------------------------------------------------
saveRDS(sce, file="mammary.rds")

## --------------------------------------------------------------------------
sessionInfo()

