## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## load vignette specific libraries
library(CRBHits)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(curl))
## compile LAST, KaKs_Calculator2.0 and DAGchainer for the vignette
vignette.paths <- make_vignette()
## set FTP handle due to travis-ci
ftp_handle <- curl::new_handle( ftp_use_epsv = FALSE,
                                crlf = TRUE,
                                ssl_verifypeer = FALSE)

## -----------------------------------------------------------------------------
## example how to check coding sequences if all are a mutiple of three

## load CDS file
cdsfile <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
cds <- Biostrings::readDNAStringSet(cdsfile)
## the following statement should return TRUE, if all sequences are a mutiple of three
all(Biostrings::width(cds) %% 3 == 0)

## -----------------------------------------------------------------------------
## example how to access CDS from URL and get longest isoform

## get coding sequences for Homo sapiens from NCBI
HOMSAP.cds.NCBI.url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/",
  "GCF_000001405.39_GRCh38.p13/",
  "GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna.gz")
HOMSAP.cds.NCBI <- Biostrings::readDNAStringSet(HOMSAP.cds.NCBI.url)
## reduce to the longest isoform
HOMSAP.cds.NCBI.longest <- isoform2longest(HOMSAP.cds.NCBI, "NCBI")
## get coding sequences for Homo sapiens from ENSEMBL
HOMSAP.cds.ENSEMBL.url <- paste0(
  "ftp://ftp.ensembl.org/pub/release-101/fasta/",
  "homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz")
HOMSAP.cds.ENSEMBL.file <- tempfile()
#download.file(HOMSAP.cds.ENSEMBL.url, HOMSAP.cds.ENSEMBL.file, quiet = FALSE)
curl::curl_download(url = HOMSAP.cds.ENSEMBL.url,
                    handle = ftp_handle,
                    mode = "wb",
                    destfile = HOMSAP.cds.ENSEMBL.file)
HOMSAP.cds.ENSEMBL <- Biostrings::readDNAStringSet(HOMSAP.cds.ENSEMBL.file)
## reduce to the longest isoform
HOMSAP.cds.ENSEMBL.longest <- isoform2longest(HOMSAP.cds.ENSEMBL, "ENSEMBL")

## get help
#?isoform2longest

## -----------------------------------------------------------------------------
## example how to get CRBHit pairs from two CDS using classical RBH

## define CDS file1
cdsfile1 <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
## define CDS file2
cdsfile2 <- system.file("fasta", "aly.cds.fasta.gz", package = "CRBHits")
## get CDS1
cds1 <- Biostrings::readDNAStringSet(cdsfile1)
## get CDS2
cds2 <- Biostrings::readDNAStringSet(cdsfile2)
## example how to perform classical RBH
ath_aly_rbh <- cds2rbh(cds1, cds2, crbh = FALSE,
                       lastpath = vignette.paths[1])
## show summary
summary(ath_aly_rbh)

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## example how to get CRBHit pairs using one thread and plot CRBHit algorithm fitting curve

## example how to perform CRBH
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE,
                        lastpath = vignette.paths[1])
## show summary
summary(ath_aly_crbh)

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## example showing cds2rbh() results

## show dimension and first parts of retained hit pairs
dim(ath_aly_crbh$crbh.pairs)
head(ath_aly_crbh$crbh.pairs)

## show first retained hit pairs for the query > target matrix
head(ath_aly_crbh$crbh1)

## get the number of CRBHit classified as rbh and sec hit pairs
table(ath_aly_crbh$crbh1$rbh_class)
table(ath_aly_crbh$crbh2$rbh_class)

## -----------------------------------------------------------------------------
## example how to use the fitting function for manual plotting

## plot fitting function
curve(ath_aly_crbh$rbh1_rbh2_fit(x),
      from = 1,
      to = 1000,
      xlab = "alnlength",
      ylab = "-log10(evalue)",
      main = "CRBH fitting")

## -----------------------------------------------------------------------------
## example how to retain single direction secondary homologs

## get CRBHit pairs with keepSingleDirection = TRUE
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE,
                        keepSingleDirection = TRUE,
                        lastpath = vignette.paths[1])
## get the number of CRBHit classified as rbh and sec hit pairs
table(ath_aly_crbh$crbh1$rbh_class)
table(ath_aly_crbh$crbh2$rbh_class)
dim(ath_aly_crbh)

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## example how to filter prior crbh for query coverage of 50%

## get CRBHit pairs with direct query coverage filter
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, qcov = 0.5,
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## detailed explanation for the Rost (1999) twilight-zone filter

## define eq2 from Rost (1999)
get_pident_by_length <- function(x){
  eq2 <- function(L){
    if(L <= 11){return(100)}
    if(L <= 450){return(480*(L^(-0.32*(1+(exp(-L/1000))))))}
    if(L > 450){return(19.5)}
  }
  return(unlist(lapply(x, eq2)))
}

## plot expected pident by alignment length using eq2 from Rost (1999)
curve(get_pident_by_length, 11, 500, pch = 20, xlab = "alignment length", ylab = "pident",
      main = "expected protein identity (eq2; Rost B. 1999)")

## -----------------------------------------------------------------------------
## example how to filter prior crbh for eq2 from Rost (1999)

## get CRBHit pairs with direct twilight-zone filter
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE, rost1999 = TRUE,
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

## get help
#?cds2rbh
## get help
#?filter.alnlen
## get help
#?filter.eval
## get help
#?filter.pident
## get help
#?filter.qcov
## get help
#?filter.rost1999
## get help
#?filter.tcov

## -----------------------------------------------------------------------------
## example for a custom filters

## define custom filter for e.g. bit score (column 12)
myfilter1 <- function(rbh, value = 500.0){
  return(dplyr::filter(rbh, bit_score >= value))
}
## define custom filter for e.g. corrected query_coverage
myfilter2 <- function(rbh, value = 0.5){
  return(dplyr::filter(rbh, ((alignment_length-mismatches-gap_opens) / query_length) >= value))
}
## get CRBHit pairs with custom filter list
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE,
                        filter = list(myfilter1, myfilter2),
                        lastpath = vignette.paths[1])
dim(ath_aly_crbh$crbh.pairs)

## -----------------------------------------------------------------------------
## example to extract CRBHit pairs classified as rbh

## reduce to rbh_class rbh
data("ath_aly_rbh", package = "CRBHits")
head(dplyr::filter(ath_aly_crbh$crbh1, rbh_class == "rbh"))
head(dplyr::filter(ath_aly_crbh$crbh2, rbh_class == "rbh"))

## -----------------------------------------------------------------------------
## example how to get crbh from two coding fasta files using median fitting

## get CRBHit pairs with median fitting
ath_aly_crbh <- cds2rbh(cds1, cds2, plotCurve = TRUE,
                        fit.type = "median",
                        lastpath = vignette.paths[1])

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## example to get a codon alignment

## define two CDS
cds1 <- Biostrings::DNAString("ATGCAACATTGC")
cds2 <- Biostrings::DNAString("ATGCATTGC")
## get codon alignment
cds2codonaln(cds1, cds2)

## get help
#?cds2codonaln

## -----------------------------------------------------------------------------
## example to alter the substitionMatrix and use the BLOSUM45 cost matrix
## for the codon alignment

## get codon alignemnt with BLOSUM45 cost matrix
cds2codonaln(cds1, cds2, substitutionMatrix = "BLOSUM45")

## -----------------------------------------------------------------------------
## example to remove codon gaps

## get codon alignment with gaps removed
cds2codonaln(cds1, cds2, remove.gaps = TRUE)

## -----------------------------------------------------------------------------
## calculate Ka/Ks on two CDS

## load example sequence data
data("ath", package="CRBHits")
data("aly", package="CRBHits")
## select a sequence pair according to a best hit pair (done for you)
cds1 <- ath[1]
cds2 <- aly[282]
## calculate Ka/Ks values on two CDS using Li model
cds2kaks(cds1, cds2, model = "Li")

## get help
#?cds2kaks

## -----------------------------------------------------------------------------
## example to use an alternative substitutionMatrix for the codon alignment
## and obtain Ka/Ks

## calculate Ka/Ks values on two CDS using Li model and BLOSUM45 cost matrix
cds2kaks(cds1, cds2, model = "Li", substitutionMatrix = "BLOSUM45")

## -----------------------------------------------------------------------------
## example how to get CRBHit pairs from two coding fasta files
cdsfile1 <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
cdsfile2 <- system.file("fasta", "aly.cds.fasta.gz", package = "CRBHits")

ath <- Biostrings::readDNAStringSet(cdsfile1)
aly <- Biostrings::readDNAStringSet(cdsfile2)

## the following function calculates CRBHit pairs using one thread and plots the fitted curve
ath_aly_crbh <- cds2rbh(cds1 = ath, cds2 = aly,
                        lastpath = vignette.paths[1])

## calculate Ka/Ks using the CRBHit pairs
ath_aly_crbh$crbh.pairs <- head(ath_aly_crbh$crbh.pairs)
ath_aly_crbh.kaks <- rbh2kaks(rbhpairs = ath_aly_crbh,
                              cds1 = ath, cds2 = aly,
                              model = "Li")
head(ath_aly_crbh.kaks)

## get help
#?rbh2kaks

## -----------------------------------------------------------------------------
## calculate Ka/Ks using the CRBHit pairs and multiple threads
ath_aly_crbh.kaks <- rbh2kaks(rbhpairs = ath_aly_crbh,
                              cds1 = ath, cds2 = aly,
                              model = "Li", threads = 2)
head(ath_aly_crbh.kaks)

## get help
#?rbh2kaks

## -----------------------------------------------------------------------------
## calculate Ka/Ks using the CRBHit pairs and multiple threads
ath_aly_crbh.kaks <- rbh2kaks(rbhpairs = ath_aly_crbh,
                              cds1 = ath, cds2 = aly,
                              model = "Li", threads = 2,
                              substitutionMatrix = "BLOSUM45")
head(ath_aly_crbh.kaks)

## get help
#?rbh2kaks

