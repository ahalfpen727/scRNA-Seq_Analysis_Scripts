## ---- echo=FALSE, results="hide", message=FALSE----------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE--------------------------------------
library(DropletUtils)

## --------------------------------------------------------------------------
# To generate the files.
example(write10xCounts, echo=FALSE) 
dir.name <- tmpdir
list.files(dir.name)

## --------------------------------------------------------------------------
sce <- read10xCounts(dir.name)
sce

## --------------------------------------------------------------------------
class(counts(sce))

## --------------------------------------------------------------------------
set.seed(1000)
mol.info.file <- DropletUtils:::sim10xMolInfo(tempfile())
mol.info.file

## --------------------------------------------------------------------------
mol.info <- read10xMolInfo(mol.info.file)
mol.info

## --------------------------------------------------------------------------
set.seed(100)
new.counts <- downsampleMatrix(counts(sce), prop=0.5)
library(Matrix)
colSums(counts(sce))
colSums(new.counts)

## --------------------------------------------------------------------------
set.seed(100)
no.sampling <- downsampleReads(mol.info.file, prop=1)
sum(no.sampling)
with.sampling <- downsampleReads(mol.info.file, prop=0.5)
sum(with.sampling)

## --------------------------------------------------------------------------
set.seed(0)
my.counts <- DropletUtils:::simCounts()

## --------------------------------------------------------------------------
br.out <- barcodeRanks(my.counts)

# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

## --------------------------------------------------------------------------
set.seed(100)
e.out <- emptyDrops(my.counts)
e.out

## --------------------------------------------------------------------------
is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)

## --------------------------------------------------------------------------
table(Limited=e.out$Limited, Significant=is.cell)

## --------------------------------------------------------------------------
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")

## --------------------------------------------------------------------------
set.seed(1000)
mult.mol.info <- DropletUtils:::sim10xMolInfo(tempfile(), nsamples=3)
mult.mol.info

## --------------------------------------------------------------------------
s.out <- swappedDrops(mult.mol.info, min.frac=0.9)
length(s.out$cleaned)
class(s.out$cleaned[[1]])

## --------------------------------------------------------------------------
sessionInfo()

