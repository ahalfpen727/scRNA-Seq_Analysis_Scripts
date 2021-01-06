## ---- echo=FALSE, results="hide"-------------------------------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir="pbmc4k")

## --------------------------------------------------------------------------
library(DropletUtils)
fname <- "pbmc4k/raw_gene_bc_matrices/GRCh38"
sce <- read10xCounts(fname, col.names=TRUE)
sce

## --------------------------------------------------------------------------
class(counts(sce))

## --------------------------------------------------------------------------
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

## --------------------------------------------------------------------------
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
    column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="MT")

## ----rankplot, fig.cap="Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts). The inferred locations of the inflection and knee points are also shown."----
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
	col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

## --------------------------------------------------------------------------
set.seed(100)
e.out <- emptyDrops(counts(sce))
sum(e.out$FDR <= 0.01, na.rm=TRUE)

## --------------------------------------------------------------------------
# using which() to automatically remove NAs.
sce <- sce[,which(e.out$FDR <= 0.01)]

## --------------------------------------------------------------------------
table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited)

## ----ambientpvalhist, fig.cap="Distribution of p-values for the assumed empty droplets."----
full.data <- read10xCounts(fname, col.names=TRUE)
set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(full.data), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
    xlab="P-value", main="", col="grey80") 

## ----qchist, fig.width=10, fig.asp=0.5, fig.cap="Histograms of QC metric distributions in the PBMC dataset."----
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
    xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
	xlab="Proportion of reads in mitochondrial genes")

## --------------------------------------------------------------------------
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
sce <- sce[,!high.mito]
summary(high.mito)

## ----abhist, fig.cap="Histogram of the log~10~-average counts for each gene in the PBMC dataset."----
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80")

## ----highexpr, fig.wide=TRUE, fig.asp=1.5, fig.cap="Percentage of total counts assigned to the top 50 most highly-abundant features in the PBMC dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell."----
plotHighestExprs(sce)

## --------------------------------------------------------------------------
library(scran)
set.seed(1000)
clusters <- quickCluster(sce, method="igraph", min.mean=0.1,
    irlba.args=list(maxit=1000)) # for convergence.
table(clusters)

## --------------------------------------------------------------------------
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))

## ----sfplot, fig.cap="Size factors for all cells in the PBMC dataset, plotted against the library size."----
plot(sce$total_counts, sizeFactors(sce), log="xy")

## --------------------------------------------------------------------------
sce <- normalize(sce)

## --------------------------------------------------------------------------
new.trend <- makeTechTrend(x=sce)

## ----trendplot, fig.cap="Variance of normalized log-expression values for each gene in the PBMC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances, while the red line represents the Poisson noise."----
fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)

## --------------------------------------------------------------------------
fit0 <- fit
fit$trend <- new.trend
dec <- decomposeVar(fit=fit)
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)

## ----hvgplot, fig.wide=TRUE, fig.cap="Distributions of normalized log-expression values for the top 10 genes with the largest biological components in the PBMC dataset. Each point represents the log-expression value in a single cell."----
plotExpression(sce, features=rownames(top.dec)[1:10])

## --------------------------------------------------------------------------
set.seed(1000)
sce <- denoisePCA(sce, technical=new.trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))

## ----screeplot, fig.cap="Variance explained by each principal component in the PBMC dataset. The red line represents the chosen number of PCs."----
plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
	ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

## ----pcaplot-init, fig.cap="Pairwise PCA plots of the first three PCs in the PBMC dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured by the log-number of expressed features.", fig.width=9----
plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")

## ----tsneplot-init, fig.cap="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to the log-number of expressed features."----
set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30)
plotTSNE(sce, colour_by="log10_total_features_by_counts")

## --------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)

## ----clustermod, fig.cap="Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links."----
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)

library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))

## ----tsneplot-cluster, fig.cap="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to its cluster identity."----
plotTSNE(sce, colour_by="Cluster")

## --------------------------------------------------------------------------
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")

## --------------------------------------------------------------------------
marker.set <- markers[["8"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

## ---- echo=FALSE, results="hide"-------------------------------------------
# Checking the cluster is what we wanted.
pf4 <- sapply(marker.set["PF4",-(1:3)], sign)
stopifnot(all(pf4==1))

## ----heatmap, fig.wide=TRUE, fig.cap="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 8 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend."----
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster),
    show_colnames=FALSE)

## --------------------------------------------------------------------------
saveRDS(sce, file="pbmc_data.rds")

## --------------------------------------------------------------------------
sessionInfo()

