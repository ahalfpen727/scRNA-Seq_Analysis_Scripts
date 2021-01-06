## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask=FALSE)
islam.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE29nnn/GSE29087/suppl/GSE29087_L139_expression_tab.txt.gz"))

## --------------------------------------------------------------------------
library(SingleCellExperiment)
counts <- read.table(islam.fname,
    colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), 
    rep("integer", 96)), skip=6, sep='\t', row.names=1)

is.spike <- grep("SPIKE", rownames(counts)) 
sce.islam <- SingleCellExperiment(list(counts=as.matrix(counts)))
isSpike(sce.islam, "spike") <- is.spike
dim(sce.islam)

## --------------------------------------------------------------------------
library(scater)
sce.islam <- calculateQCMetrics(sce.islam)
sce.islam$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))

libsize.drop <- isOutlier(sce.islam$total_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce.islam$grouping)
feature.drop <- isOutlier(sce.islam$total_features_by_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce.islam$grouping)
spike.drop <- isOutlier(sce.islam$pct_counts_spike, nmads=3, type="higher", 
    batch=sce.islam$grouping)
    
sce.islam <- sce.islam[,!(libsize.drop | feature.drop | 
    spike.drop | sce.islam$grouping=="Neg")]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce.islam))

## --------------------------------------------------------------------------
library(scran)
sce.islam <- computeSpikeFactors(sce.islam, general.use=TRUE)
head(sizeFactors(sce.islam))    
head(sizeFactors(sce.islam, "spike")) # same as general size factors.

## --------------------------------------------------------------------------
sce.islam <- normalize(sce.islam)

## --------------------------------------------------------------------------
deconv.sf <- computeSumFactors(sce.islam, sf.out=TRUE, cluster=sce.islam$grouping)
head(deconv.sf)

## ----normplotspikemef, fig.cap="Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the mESC/MEF dataset. Axes are shown on a log-scale, and cells are coloured according to their identity. Deconvolution size factors were computed with small pool sizes owing to the low number of cells of each type."----
colours <- c(mESC="red", MEF="grey")
plot(sizeFactors(sce.islam), deconv.sf, col=colours[sce.islam$grouping], pch=16, 
    log="xy", xlab="Size factor (spike-in)", ylab="Size factor (deconvolution)")
legend("bottomleft", col=colours, legend=names(colours), pch=16)

## --------------------------------------------------------------------------
sessionInfo()

