## ---- echo=FALSE, results="hide", message=FALSE----------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----style, echo=FALSE, results='asis'-------------------------------------
BiocStyle::markdown()

## --------------------------------------------------------------------------
library(TENxBrainData)
tenx <- TENxBrainData()
tenx

## --------------------------------------------------------------------------
counts(tenx)

## --------------------------------------------------------------------------
options(DelayedArray.block.size=2e9)

## ---- eval = FALSE---------------------------------------------------------
#  lib.sizes <- colSums(counts(tenx))
#  n.exprs <- colSums(counts(tenx) != 0L)
#  ave.exprs <- rowMeans(counts(tenx))

## --------------------------------------------------------------------------
tenx20k <- tenx[, seq_len(20000)]
chunksize <- 5000
cidx <- snow::splitIndices(ncol(tenx20k), ncol(tenx20k) / chunksize)

## --------------------------------------------------------------------------
lib.sizes <- n.exprs <- numeric(ncol(tenx20k))
tot.exprs <- numeric(nrow(tenx20k))
for (i in head(cidx, 2)) {
    message(".", appendLF=FALSE)
    m <- as.matrix(counts(tenx20k)[,i, drop=FALSE])
    lib.sizes[i] <- colSums(m)
    n.exprs[i] <- colSums(m != 0)
    tot.exprs <- tot.exprs + rowSums(m)
    }
ave.exprs <- tot.exprs / ncol(tenx20k)

## --------------------------------------------------------------------------
colData(tenx20k)$lib.sizes <- lib.sizes
colData(tenx20k)$n.exprs <- n.exprs
rowData(tenx20k)$ave.exprs <- ave.exprs

## --------------------------------------------------------------------------
hist(
    log10(colData(tenx20k)$lib.sizes),
    xlab=expression(Log[10] ~ "Library size"),
    col="grey80"
)

## --------------------------------------------------------------------------
hist(colData(tenx20k)$n.exprs, xlab="Number of detected genes", col="grey80")

## --------------------------------------------------------------------------
hist(
    log10(rowData(tenx20k)$ave.exprs),
    xlab=expression(Log[10] ~ "Average count"),
    col="grey80"
)

## --------------------------------------------------------------------------
o <- order(rowData(tenx20k)$ave.exprs, decreasing=TRUE)
head(rowData(tenx20k)[o,])

## ---- eval=FALSE-----------------------------------------------------------
#  destination <- tempfile()
#  saveRDS(tenx, file = destination)

## --------------------------------------------------------------------------
library(BiocParallel)
register(bpstart(SnowParam(5)))

## --------------------------------------------------------------------------
iterator <- function(tenx, cols_per_chunk = 5000, n = Inf) {
    start <- seq(1, ncol(tenx), by = cols_per_chunk)
    end <- c(tail(start, -1) - 1L, ncol(tenx))
    n <- min(n, length(start))
    i <- 0L
    function() {
        if (i == n)
            return(NULL)
        i <<- i + 1L
        c(start[i], end[i])
    }
}

## --------------------------------------------------------------------------
iter <- iterator(tenx)
iter()
iter()
iter()

## --------------------------------------------------------------------------
fun <- function(crng, counts, ...) {
    ## `fun()` needs to be self-contained for some parallel back-ends
    suppressPackageStartupMessages({
        library(TENxBrainData)
    })
    m <- as.matrix( counts[ , seq(crng[1], crng[2]) ] )
    list(
        row = list(
            n = rowSums(m != 0), sum = rowSums(m), sumsq = rowSums(m * m)
        ),
        column = list(
            n = colSums(m != 0), sum = colSums(m), sumsq = colSums(m * m)
        )
    )
}

## --------------------------------------------------------------------------
res <- fun( iter(), unname(counts(tenx)) )
str(res)

## --------------------------------------------------------------------------
reduce <- function(x, y) {
    list(
        row = Map(`+`, x$row, y$row),
        column = Map(`c`, x$column, y$column)
    )
}

## --------------------------------------------------------------------------
str( reduce(res, res) )

## --------------------------------------------------------------------------
res <- bpiterate(
    iterator(tenx, n = 5), fun, counts = unname(counts(tenx)), 
    REDUCE = reduce, reduce.in.order = TRUE
)
str(res)

## --------------------------------------------------------------------------
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## --------------------------------------------------------------------------
h5ls(fname)

## --------------------------------------------------------------------------
start <- h5read(fname, "/mm10/indptr", start=1, count=25001)
head(start)

## --------------------------------------------------------------------------
library(data.table)
dt <- data.table(
    row = h5read(fname, "/mm10/indices", start = 1, count = tail(start, 1)) + 1,
    column = rep(seq_len(length(start) - 1), diff(start)),
    count = h5read(fname, "/mm10/data", start = 1, count = tail(start, 1))
)
dt

## --------------------------------------------------------------------------
dt[ , 
    list(n = .N, sum = sum(count), sumsq = sum(count * count)),
    keyby=row]
dt[ , 
    list(n = .N, sum = sum(count), sumsq = sum(count * count)),
    keyby=column]

## --------------------------------------------------------------------------
sessionInfo()

