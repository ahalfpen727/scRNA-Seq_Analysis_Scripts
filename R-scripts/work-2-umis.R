## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
base.url <- file.path("https://storage.googleapis.com",
    "linnarsson-lab-www-blobs/blobs/cortex")
mRNA.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_mRNA_17-Aug-2014.txt"))
mito.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_mito_17-Aug-2014.txt"))
spike.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_spikes_17-Aug-2014.txt"))

## --------------------------------------------------------------------------
readFormat <- function(infile) { 
    # First column is empty.
    metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    metadata <- as.data.frame(t(metadata))

    # First column after row names is some useless filler.
    counts <- read.delim(infile, stringsAsFactors=FALSE, 
        header=FALSE, row.names=1, skip=11)[,-1] 
    counts <- as.matrix(counts)
    return(list(metadata=metadata, counts=counts))
}

## --------------------------------------------------------------------------
endo.data <- readFormat(mRNA.path)
spike.data <- readFormat(spike.path)
mito.data <- readFormat(mito.path)

## --------------------------------------------------------------------------
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]

## ---- echo=FALSE-----------------------------------------------------------
stopifnot(identical(endo.data$metadata$cell_id, spike.data$metadata$cell_id)) # should be the same.
stopifnot(all(endo.data$metadata$cell_id==mito.data$metadata$cell_id)) # should now be the same.

## --------------------------------------------------------------------------
raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
endo.data$counts <- new.counts

## --------------------------------------------------------------------------
library(SingleCellExperiment)
all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)
dim(sce)

## --------------------------------------------------------------------------
# Specifying the nature of each row.
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
isSpike(sce, "Spike") <- is.spike

# Adding Ensembl IDs.
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
rowData(sce)$ENSEMBL <- ensembl

sce

## ---- echo=FALSE, results='hide'-------------------------------------------
# Save some memory.
rm(mito.data, endo.data, spike.data, new.counts)
gc()

## --------------------------------------------------------------------------
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

## ----libplotbrain, fig.wide=TRUE, fig.cap="Histograms of QC metrics including the library sizes, number of expressed genes and proportion of UMIs assigned to spike-in transcripts or mitochondrial genes for all cells in the brain dataset."----
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)",
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")

## --------------------------------------------------------------------------
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")

## --------------------------------------------------------------------------
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
    BySpike=sum(spike.drop), Remaining=ncol(sce))

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## ----phaseplotbrain, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the brain dataset, where each point represents a cell."----
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
table(assignments$phase)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## ----topgenebrain, fig.asp=1.2, fig.wide=TRUE, fig.cap="Percentage of total counts assigned to the top 50 most highly-abundant features in the brain dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature."----
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize

## ----abhistbrain, fig.cap="Histogram of log-average counts for all genes in the brain dataset."----
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey",
    xlab=expression(Log[10]~"average count"))

## --------------------------------------------------------------------------
rowData(sce)$ave.count <- ave.counts
to.keep <- ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## --------------------------------------------------------------------------
set.seed(1000)
clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
summary(sizeFactors(sce))

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## ----normplotbrain, fig.cap="Size factors from deconvolution, plotted against library sizes for all cells in the brain dataset. Axes are shown on a log-scale."----
plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
    ylab="Library size (thousands)", xlab="Size factor")

## --------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, type="Spike", general.use=FALSE)

## --------------------------------------------------------------------------
sce <- normalize(sce)

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## --------------------------------------------------------------------------
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.4))
var.out <- decomposeVar(sce, var.fit)

## ----hvgplotbrain, fig.cap="Variance of normalized log-expression values against the mean for each gene, calculated across all cells in the brain dataset after blocking on the sex effect. The blue line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points)."----
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
points(var.out$mean[isSpike(sce)], var.out$total[isSpike(sce)], col="red", pch=16)
curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

## ----hvgvioplotbrain, fig.cap="Violin plots of normalized log-expression values for the top 10 HVGs in the brain dataset. For each gene, each point represents the log-expression value for an individual cell."----
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes], 
    point_alpha=0.05, jitter_type="jitter") + fontsize

## --------------------------------------------------------------------------
set.seed(1000)
sce <- denoisePCA(sce, technical=var.fit$trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))

## ---- echo=FALSE, results='hide', message=FALSE----------------------------
gc()

## ----tsneplotbrain, fig.cap="_t_-SNE plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of _Neurod6_ (left) or _Mog_ (right).", fig.width=12, fig.asp=0.5----
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50)
tsne1 <- plotTSNE(sce, colour_by="Neurod6") + fontsize
tsne2 <- plotTSNE(sce, colour_by="Mog") + fontsize
multiplot(tsne1, tsne2, cols=2)

## ----pcaplotbrain, fig.cap="PCA plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of the _Neurod6_ (left) or _Mog_ (right).", fig.width=12, fig.asp=0.5----
pca1 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Neurod6") + fontsize
pca2 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Mog") + fontsize
multiplot(pca1, pca2, cols=2)

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## --------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters <- cluster.out$membership
table(my.clusters)

## ----tsneclusterbrain, message=FALSE, fig.cap="_t_-SNE plot of the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity."----
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## ----fdlbrain, message=FALSE, fig.cap="Force-directed layout for the shared nearest-neighbour graph of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity."----
set.seed(2000)
reducedDim(sce, "force") <- igraph::layout_with_fr(snn.gr, niter=5000)
plotReducedDim(sce, colour_by="cluster", use_dimred="force")

## --------------------------------------------------------------------------
igraph::modularity(cluster.out)

## ----heatmodbrain, fig.cap="Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links."----
mod.out <- clusterModularity(snn.gr, my.clusters, get.values=TRUE)
ratio <- mod.out$observed/mod.out$expected
lratio <- log10(ratio + 1)

library(pheatmap)
pheatmap(lratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))

## ----graphbrain, fig.cap="Force-directed layout showing the relationships between clusters based on the ratio of observed to expected total weights between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding ratio."----
cluster.gr <- igraph::graph_from_adjacency_matrix(ratio, 
    mode="undirected", weighted=TRUE, diag=FALSE)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*10)  

## ---- echo=FALSE, results="hide"-------------------------------------------
old.digits <- options()$digits
options(digits=3)

## --------------------------------------------------------------------------
markers <- findMarkers(sce, my.clusters, direction="up")
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

## ---- echo=FALSE, results="hide"-------------------------------------------
# Checking the cluster is what we wanted, along with cluster 10 (=9 in marker.set).
gad1 <- sapply(marker.set["Gad1",-(1:3)], sign)
stopifnot(gad1[9]==-1)
stopifnot(all(gad1[-9]==1))

gad2 <- sapply(marker.set["Gad2",-(1:3)], sign)
stopifnot(gad2[9]==-1)
stopifnot(all(gad2[-9]==1))

stopifnot(all(sapply(marker.set["Synpr",-(1:3)], sign)==1))

options(digits=old.digits)

## --------------------------------------------------------------------------
gzout <- gzfile("brain_marker_1.tsv.gz", open="wb")
write.table(marker.set, file=gzout, sep="\t", quote=FALSE, col.names=NA)
close(gzout)

## ----heatmapmarkerbrain, fig.wide=TRUE, fig.cap="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the brain dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend."----
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(my.clusters),
    colour_columns_by="cluster", cluster_cols=FALSE, 
    center=TRUE, symmetric=TRUE, zlim=c(-5, 5))

## --------------------------------------------------------------------------
saveRDS(file="brain_data.rds", sce)

## ---- echo=FALSE, results='hide'-------------------------------------------
gc()

## --------------------------------------------------------------------------
sessionInfo()

