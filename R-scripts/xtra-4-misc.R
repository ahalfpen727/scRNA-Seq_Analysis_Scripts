## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
library(scran)
sce.hsc <- readRDS("hsc_data.rds")    

set.seed(100)
var.cor <- correlatePairs(sce.hsc, subset.row=grep("^H2-", rownames(sce.hsc)))
head(var.cor)

## --------------------------------------------------------------------------
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)

## --------------------------------------------------------------------------
correlatePairs(sce.hsc, subset.row=cbind("Fos", "Jun"))

## ----fosjuncorplot, fig.cap="Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC dataset."----
library(scater)
plotExpression(sce.hsc, features="Fos", x="Jun")

## --------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
mahata.fname <- bfcrpath(bfc, 
    "http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx")

## --------------------------------------------------------------------------
library(readxl)
incoming <- as.data.frame(read_excel(mahata.fname, sheet=1))
rownames(incoming) <- incoming[,1]
incoming <- incoming[,-1]
incoming <- incoming[,!duplicated(colnames(incoming))] # Remove duplicated genes.
sce.th2 <- SingleCellExperiment(list(logcounts=t(incoming)))

## ----phaseplotth2, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell."----
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce.th2), keytype="SYMBOL", column="ENSEMBL")

set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))
assignments <- cyclone(sce.th2, mm.pairs, gene.names=ensembl, assay.type="logcounts")

plot(assignments$score$G1, assignments$score$G2M, 
    xlab="G1 score", ylab="G2/M score", pch=16)

## --------------------------------------------------------------------------
design <- model.matrix(~ G1 + G2M, assignments$score)
fit.block <- trendVar(sce.th2, design=design, parametric=TRUE, use.spikes=NA)
dec.block <- decomposeVar(sce.th2, fit.block)

library(limma)
sce.th2.block <- sce.th2
assay(sce.th2.block, "corrected") <- removeBatchEffect(
    logcounts(sce.th2), covariates=design[,-1])

sce.th2.block <- denoisePCA(sce.th2.block, technical=dec.block, 
    assay.type="corrected")
dim(reducedDim(sce.th2.block, "PCA"))

## ----pcaplotth2, fig.width=12, fig.asp=0.5, fig.cap="PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively."----
sce.th2$G1score <- sce.th2.block$G1score <- assignments$score$G1
sce.th2$G2Mscore <- sce.th2.block$G2Mscore <- assignments$score$G2M

# Without blocking on phase score.
fit <- trendVar(sce.th2, parametric=TRUE, use.spikes=NA) 
sce.th2 <- denoisePCA(sce.th2, technical=fit$trend)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
out <- plotReducedDim(sce.th2, use_dimred="PCA", ncomponents=2, colour_by="G1score", 
    size_by="G2Mscore") + fontsize + ggtitle("Before removal")

# After blocking on the phase score.
out2 <- plotReducedDim(sce.th2.block, use_dimred="PCA", ncomponents=2, 
    colour_by="G1score", size_by="G2Mscore") + fontsize + 
    ggtitle("After removal")
multiplot(out, out2, cols=2)

## ----diffusionth2, fig.cap="A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_. A larger `sigma` is used compared to the default value to obtain a smoother plot."----
plotDiffusionMap(sce.th2.block, colour_by="Gata3",
    run_args=list(use_dimred="PCA", sigma=25)) + fontsize

## --------------------------------------------------------------------------
sessionInfo()

