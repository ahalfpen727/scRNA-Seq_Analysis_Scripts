## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(TENxBrainData)
sce <- TENxBrainData20k() # downloads once and caches it for future use.
sce

## --------------------------------------------------------------------------
counts(sce)
object.size(counts(sce))
file.info(path(counts(sce)))$size

## --------------------------------------------------------------------------
tmp <- counts(sce)
tmp <- log2(tmp + 1)
tmp

## --------------------------------------------------------------------------
library(scater)
sce <- calculateQCMetrics(sce, compact=TRUE) # compacting for clean output.
sce$scater_qc

## --------------------------------------------------------------------------
bpp <- MulticoreParam(2)
bpp

## --------------------------------------------------------------------------
bpp <- SnowParam(5)
bpp

## ---- eval=FALSE-----------------------------------------------------------
#  bpp <- BatchtoolsParam(10, cluster="slurm",
#  	resources=list(walltime=20000, memory=8000, ncpus=1))

## --------------------------------------------------------------------------
alt <- calculateQCMetrics(sce, BPPARAM=MulticoreParam(2), compact=TRUE)

## ---- echo=FALSE-----------------------------------------------------------
if (!isTRUE(all.equal(alt, sce))) {
	stop("parallelization changes the result")
}

## --------------------------------------------------------------------------
all.equal(alt, sce) 

## --------------------------------------------------------------------------
sce.pbmc <- readRDS("pbmc_data.rds")

## --------------------------------------------------------------------------
library(scran)
library(BiocNeighbors)
snn.gr <- buildSNNGraph(sce.pbmc, BNPARAM=AnnoyParam(), use.dimred="PCA")

## --------------------------------------------------------------------------
clusters <- igraph::cluster_walktrap(snn.gr)
table(Exact=sce.pbmc$Cluster, Approx=clusters$membership)

## --------------------------------------------------------------------------
sessionInfo()

