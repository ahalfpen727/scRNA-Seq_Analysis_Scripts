## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
lun.zip <- bfcrpath(bfc, 
    file.path("https://www.ebi.ac.uk/arrayexpress/files",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"))
lun.sdrf <- bfcrpath(bfc, 
    file.path("https://www.ebi.ac.uk/arrayexpress/files",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))
unzip(lun.zip)

## --------------------------------------------------------------------------
plate1 <- read.delim("counts_Calero_20160113.tsv", 
    header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim("counts_Calero_20160325.tsv", 
    header=TRUE, row.names=1, check.names=FALSE)

gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))

## --------------------------------------------------------------------------
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

## --------------------------------------------------------------------------
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce

## --------------------------------------------------------------------------
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))

## --------------------------------------------------------------------------
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)

## --------------------------------------------------------------------------
metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))

## --------------------------------------------------------------------------
colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)

## --------------------------------------------------------------------------
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

## --------------------------------------------------------------------------
new.names <- rowData(sce)$SYMBOL
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ENSEMBL[missing.name]
dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
head(rownames(sce))

## --------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

## --------------------------------------------------------------------------
library(scater)
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)

## ----qcplot416b, fig.wide=TRUE, fig.cap="Distributions of various QC metrics for all cells in the 416B dataset. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes."----
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(
    plotColData(sce, y="total_counts", x="PlateOnco"),
    plotColData(sce, y="total_features_by_counts", x="PlateOnco"),
    plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
    plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
    cols=2)

## ----qcbiplot416b, fig.width=10, fig.asp=0.5, fig.cap="Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the 416B dataset."----
par(mfrow=c(1,3))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes",
    ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_ERCC, xlab="Number of expressed genes",
    ylab="ERCC proportion (%)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes",
    ylab="Mitochondrial proportion (%)")

## --------------------------------------------------------------------------
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$PlateOnco)

## --------------------------------------------------------------------------
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
    batch=sce$PlateOnco)

## --------------------------------------------------------------------------
keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=sum(keep))

## --------------------------------------------------------------------------
sce$PassQC <- keep
saveRDS(sce, file="416B_preQC.rds")
sce <- sce[,keep]
dim(sce)

## --------------------------------------------------------------------------
attr(libsize.drop, "thresholds")
attr(spike.drop, "thresholds")

## --------------------------------------------------------------------------
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)

## ----phaseplot416b, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases."----
plot(assignments$score$G1, assignments$score$G2M, 
    xlab="G1 score", ylab="G2/M score", pch=16)

## --------------------------------------------------------------------------
sce$phases <- assignments$phases
table(sce$phases)

## ----topgene416b, fig.asp=1.2, fig.wide=TRUE, fig.cap="Percentage of total counts assigned to the top 50 most highly-abundant features in the 416B dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature."----
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize

## ----abhist416b, fig.cap="Histogram of log-average counts for all genes in the 416B dataset."----
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))

## --------------------------------------------------------------------------
demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)

## ----nexprshist416b, fig.cap="The number of cells expressing each gene in the 416B dataset, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location."----
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
    xlab=expression(Log[10]~"average count"))

## --------------------------------------------------------------------------
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)

## --------------------------------------------------------------------------
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

## ----normplot416b, fig.cap="Size factors from deconvolution, plotted against library sizes for all cells in the 416B dataset. Axes are shown on a log-scale. Wild-type cells are shown in black and oncogene-induced cells are shown in red."----
plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
    xlab="Library size (millions)", ylab="Size factor",
    col=c("red", "black")[sce$Oncogene], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
    legend=levels(sce$Oncogene))

## --------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)

## --------------------------------------------------------------------------
sce <- normalize(sce)

## ---- echo=FALSE, results="hide"-------------------------------------------
gc()

## --------------------------------------------------------------------------
var.fit <- trendVar(sce, parametric=TRUE, block=sce$Plate,
    loess.args=list(span=0.3))
var.out <- decomposeVar(sce, var.fit)
head(var.out)

## ----hvgplot416b, fig.cap="Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

## ----hvgvioplot416b, fig.cap="Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the 416B dataset. Each point represents the log-expression value in a single cell."----
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize

## --------------------------------------------------------------------------
library(limma)
assay(sce, "corrected") <- removeBatchEffect(logcounts(sce), 
    design=model.matrix(~sce$Oncogene), batch=sce$Plate)
assayNames(sce)

## --------------------------------------------------------------------------
sce <- denoisePCA(sce, technical=var.out, assay.type="corrected")
dim(reducedDim(sce, "PCA")) 

## --------------------------------------------------------------------------
sce2 <- denoisePCA(sce, technical=var.fit$trend, 
    assay.type="corrected", value="lowrank") 
assayNames(sce2)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(sce2)
gc()

## ----pcaplot416b-onco, fig.cap="Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to oncogene induction status.", fig.width=9----
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
    colour_by="Oncogene") + fontsize

## ----pcaplot416b-batch, fig.cap="Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to the plate of origin.", fig.width=9----
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
    colour_by="Plate") + fontsize

## ----tsneplot416b, fig.cap="_t_-SNE plots constructed from the denoised PCs in the 416B dataset, using a range of perplexity values. Each point represents a cell, coloured according to its oncogene induction status. Bars represent the coordinates of the cells on each axis.", fig.width=15, fig.asp=0.5----
set.seed(100)
out5 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=5),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 5")

set.seed(100)
out10 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=10),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 10")

set.seed(100)
out20 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=20),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 20")

multiplot(out5, out10, out20, cols=3)

## --------------------------------------------------------------------------
set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=20)
reducedDimNames(sce)

## --------------------------------------------------------------------------
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")

## --------------------------------------------------------------------------
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
    minClusterSize=10, verbose=0))

## --------------------------------------------------------------------------
table(my.clusters, sce$Plate)
table(my.clusters, sce$Oncogene)

## ----tsnecluster416b, fig.cap="_t_-SNE plot of the denoised PCs of the 416B dataset. Each point represents a cell and is coloured according to the cluster identity to which it was assigned."----
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize

## ----silhouette416b, fig.cap="Barplot of silhouette widths for cells in each cluster. Each cluster is assigned a colour and cells with positive widths are coloured according to the colour of its assigned cluster. Any cell with a negative width is coloured according to the colour of the cluster that it is closest to. The average width for all cells in each cluster is shown, along with the average width for all cells in the dataset."----
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
    border=sil.cols, col=sil.cols, do.col.sort=FALSE) 

## ----echo=FALSE, results='hide'--------------------------------------------
gc()

## --------------------------------------------------------------------------
markers <- findMarkers(sce, my.clusters, block=sce$Plate)

## ---- echo=FALSE, results="hide"-------------------------------------------
old.digits <- options()$digits
options(digits=3)

## --------------------------------------------------------------------------
marker.set <- markers[["1"]]
head(marker.set, 10)

## ---- echo=FALSE, results="hide"-------------------------------------------
# Crashing if cluster 1 is not what we think it is; 
# thus, avoid mismatch between text and results.
stopifnot(all(sapply(marker.set["Myh11",-(1:3)], sign)==1))
stopifnot(all(sapply(marker.set["Mcm2",-(1:3)], sign)==-1))

options(digits=old.digits)

## --------------------------------------------------------------------------
write.table(marker.set, file="416B_marker_1.tsv", sep="\t", 
    quote=FALSE, col.names=NA)

## ----heatmapmarker416b, fig.width=10, fig.asp=0.8, fig.cap="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, the plate of origin or the oncogene induction status of each cell, as indicated by the legend."----
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
    colour_columns_by=c("cluster", "Plate", "Oncogene"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5)) 

## --------------------------------------------------------------------------
saveRDS(file="416B_data.rds", sce)

## ---- echo=FALSE, results='hide'-------------------------------------------
gc()

## --------------------------------------------------------------------------
sessionInfo()

