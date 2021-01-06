## ----style, echo=FALSE, results='asis'-------------------------------------
BiocStyle::markdown()

## ---- message=FALSE--------------------------------------------------------
library(scRNAseq)
data(fluidigm)
fluidigm

## --------------------------------------------------------------------------
head(assay(fluidigm)[,1:3])

## --------------------------------------------------------------------------
head(assay(fluidigm, 2)[,1:3])

## --------------------------------------------------------------------------
names(assays(fluidigm))
head(assays(fluidigm)$fpkm[,1:3])

## --------------------------------------------------------------------------
dim(fluidigm)
table(rowSums(assay(fluidigm))>0)

## --------------------------------------------------------------------------
colData(fluidigm)

## --------------------------------------------------------------------------
names(metadata(fluidigm))
metadata(fluidigm)$which_qc

## --------------------------------------------------------------------------
data(th2)
ercc_idx <- grep("^ERCC-", rownames(th2))
th2_endogenous <- th2[-ercc_idx,]
th2_ercc <- th2[ercc_idx,]

head(assay(th2_ercc)[,1:4])

