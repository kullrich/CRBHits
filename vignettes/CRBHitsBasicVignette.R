## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(CRBHits)
#compile LAST and KaKs_Calculator2.0 for the vignette
vignette.paths <- make.vignette()
#example how to check coding sequences if all are a mutiple of three
#cdsfile <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
#cds <- Biostrings::readDNAStringSet(cdsfile)
data("ath")
cds <- ath
#the following statement should return TRUE, if all sequences are a mutiple of three
all(Biostrings::width(cds) %% 3 == 0)

## -----------------------------------------------------------------------------
#example how to get crbh from two coding fasta files
#cdsfile1 <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
#cdsfile2 <- system.file("fasta", "aly.cds.fasta.gz", package = "CRBHits")
#cds1 <- Biostrings::readDNAStringSet(cdsfile1)
#cds2 <- Biostrings::readDNAStringSet(cdsfile2)
data("ath")
data("aly")
cds1 <- ath
cds2 <- aly
#the following function calculates crbh matrix using one thread and plots the fitted curve
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE,
                        lastpath = vignette.paths[1])
summary(ath_aly_crbh)

## -----------------------------------------------------------------------------
#example how to perform classical rbh
ath_aly_rbh <- cds2rbh(cds1, cds2, crbh = FALSE,
                       lastpath = vignette.paths[1])

## -----------------------------------------------------------------------------
#show dimension and first retained hit pairs
dim(ath_aly_crbh$crbh.pairs)
head(ath_aly_crbh$crbh.pairs)
#show first retained hit pairs for the query > target matrix
head(ath_aly_crbh$crbh1)
#get the number of rbh and sec hit pairs
table(ath_aly_crbh$crbh1$rbh_class)
table(ath_aly_crbh$crbh2$rbh_class)

## -----------------------------------------------------------------------------
#example how to use the fitting function for manual plotting
curve(ath_aly_crbh$rbh1_rbh2_fit(x),
      from = 1,
      to = 1000,
      xlab = "alnlength",
      ylab = "-log10(evalue)",
      main = "CRBH fitting")

## -----------------------------------------------------------------------------
#example how to retain single direction secondary homologs
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, keepSingleDirection = TRUE,
                        lastpath = vignette.paths[1])
table(ath_aly_crbh$crbh1$rbh_class)
table(ath_aly_crbh$crbh2$rbh_class)
dim(ath_aly_crbh)

## -----------------------------------------------------------------------------
#example how to filter prior crbh for query coverage
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, qcov = 0.5,
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

## -----------------------------------------------------------------------------
#plot expected pident by alignment length using eq2 from Rost (1999)
get_pident_by_length <- function(x){
  eq2 <- function(L){
    if(L <= 11){return(100)}
    if(L <= 450){return(480*(L^(-0.32*(1+(exp(-L/1000))))))}
    if(L > 450){return(19.5)}
  }
  return(unlist(lapply(x, eq2)))
}
curve(get_pident_by_length, 11, 500, pch = 20, xlab = "alignment length", ylab = "pident",
      main = "expected protein identity (eq2; Rost B. 1999)")

## -----------------------------------------------------------------------------
#example how to filter prior crbh for eq2 from Rost (1999)
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, rost1999 = TRUE,
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

## -----------------------------------------------------------------------------
#example for a custom filter for e.g. bit score (column 12)
myfilter <- function(rbh, value = 500.0){
  return(rbh[as.numeric(rbh[, 12]) >= value , , drop = FALSE])
}
#example hot to filter prior crbh with custom filter
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, filter = list(myfilter),
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

