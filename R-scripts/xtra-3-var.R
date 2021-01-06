## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask=FALSE)
wilson.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE61nnn/GSE61533/suppl/GSE61533_HTSEQ_count_results.xls.gz"))

library(R.utils)
wilson.name2 <- "GSE61533_HTSEQ_count_results.xls"
gunzip(wilson.fname, destname=wilson.name2, remove=FALSE, overwrite=TRUE)

## --------------------------------------------------------------------------
library(gdata)
all.counts <- read.xls(wilson.name2, sheet=1, header=TRUE)
rownames(all.counts) <- all.counts$ID
all.counts <- as.matrix(all.counts[,-1])

## --------------------------------------------------------------------------
library(SingleCellExperiment)
sce.hsc <- SingleCellExperiment(list(counts=all.counts))
dim(sce.hsc)
is.spike <- grepl("^ERCC", rownames(sce.hsc))
isSpike(sce.hsc, "ERCC") <- is.spike
summary(is.spike)

## --------------------------------------------------------------------------
library(scater)
sce.hsc <- calculateQCMetrics(sce.hsc)
libsize.drop <- isOutlier(sce.hsc$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce.hsc$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce.hsc$pct_counts_ERCC, nmads=3, type="higher")
sce.hsc <- sce.hsc[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce.hsc))

## --------------------------------------------------------------------------
to.keep <- nexprs(sce.hsc, byrow=TRUE) > 0
sce.hsc <- sce.hsc[to.keep,]
summary(to.keep)

## ---- warning=FALSE--------------------------------------------------------
library(scran)
sce.hsc <- computeSumFactors(sce.hsc)
summary(sizeFactors(sce.hsc))
sce.hsc <- computeSpikeFactors(sce.hsc, type="ERCC", general.use=FALSE)
summary(sizeFactors(sce.hsc, "ERCC"))
sce.hsc <- normalize(sce.hsc)

## ----hvgplothsc, fig.cap="Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
var.fit <- trendVar(sce.hsc, parametric=TRUE, loess.args=list(span=0.3))
var.out <- decomposeVar(sce.hsc, var.fit)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce.hsc)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

## --------------------------------------------------------------------------
hvg.out <- var.out[which(var.out$FDR <= 0.05),]
nrow(hvg.out)

## --------------------------------------------------------------------------
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
head(hvg.out)

## ----hvgvioplothsc, fig.cap="Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the HSC dataset. Each point represents the log-expression value in a single cell."----
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotExpression(sce.hsc, features=rownames(hvg.out)[1:10]) + fontsize

## --------------------------------------------------------------------------
saveRDS(sce.hsc, file="hsc_data.rds")

## --------------------------------------------------------------------------
var.fit.nospike <- trendVar(sce.hsc, parametric=TRUE, 
    use.spikes=FALSE, loess.args=list(span=0.2))
var.out.nospike <- decomposeVar(sce.hsc, var.fit.nospike)

## ----hvgplot416b2, fig.cap="Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes (black), with spike-in transcripts shown in red."----
plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)

## ----trendplotblock-416b, fig.cap="Plate-specific variance estimates for all spike-in transcripts in the 416B dataset, plotted against the plate-specific means. Each point represents a spike-in transcript, numbered by the plate from which the values were estimated. The red line denotes the fitted mean-variance trend."----
# Loading the saved object.
sce.416B <- readRDS("416B_data.rds") 

# Repeating the trendVar() call.
var.fit <- trendVar(sce.416B, parametric=TRUE, block=sce.416B$Plate,
    loess.args=list(span=0.3))

matplot(var.fit$means, var.fit$vars, col=c("darkorange", "forestgreen"),
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit$trend(x), add=TRUE, col="red")

## ----sizefacplot-416b, fig.width=10, fig.asp=0.5, fig.cap="Plate-specific distribution of the size factors for endogenous genes (left) and spike-in transcripts (right)."----
tmp.416B <- sce.416B
tmp.416B$log_size_factor <- log(sizeFactors(sce.416B))
tmp.416B$log_size_factor_ERCC <- log(sizeFactors(sce.416B, "ERCC"))
p1 <- plotColData(tmp.416B, x="Plate", y="log_size_factor")
p2 <- plotColData(tmp.416B, x="Plate", y="log_size_factor_ERCC")
multiplot(p1, p2, cols=2)

## --------------------------------------------------------------------------
sce.416B.2 <- multiBlockNorm(sce.416B, sce.416B$Plate)
comb.out <- multiBlockVar(sce.416B.2, block=sce.416B.2$Plate,
    trend.args=list(parametric=TRUE, loess.args=list(span=0.4)))

## --------------------------------------------------------------------------
head(comb.out[,1:6])

## ----hvgplotbatch416b, fig.width=10, fig.asp=0.5, fig.cap="Variance of normalized log-expression values for each gene in each plate of the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
par(mfrow=c(1,2))
is.spike <- isSpike(sce.416B.2)
for (plate in levels(sce.416B.2$Plate)) {
    cur.out <- comb.out$per.block[[plate]]
    plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
        ylab="Variance of log-expression", main=plate)
    curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}

## --------------------------------------------------------------------------
lfit <- trendVar(sce.416B, design=model.matrix(~sce.416B$Plate))

## --------------------------------------------------------------------------
sessionInfo()

