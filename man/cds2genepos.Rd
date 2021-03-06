% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cds2genepos.R
\name{cds2genepos}
\alias{cds2genepos}
\title{cds2genepos}
\usage{
cds2genepos(cds, source = "NCBI", keep.names = NULL)
}
\arguments{
\item{cds}{\code{DNAStringSet} [mandatory]}

\item{source}{source indicating either NCBI or ENSEMBL [default: NCBI]}

\item{keep.names}{vector indicating gene ids to be kept before chromosomal position assignment [default: NULL]}
}
\value{
\code{matrix}
1: $gene.seq.id\cr
2: $gene.chr\cr
3: $gene.start\cr
4: $gene.end\cr
5: $gene.mid\cr
6: $gene.strand\cr
7: $gene.idx\cr
}
\description{
This function extracts the gene position from either NCBI or ENSEMBL CDS input.
}
\examples{
\dontrun{
## load example sequence data
## set EnsemblPlants URL
ensemblPlants <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/"
## set Arabidopsis thaliana CDS URL
ARATHA.cds.url <- paste0(ensemblPlants,
 "arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz")
ARATHA.cds.file <- tempfile()
## download CDS
download.file(ARATHA.cds.url, ARATHA.cds.file, quiet = FALSE)
ARATHA.cds <- Biostrings::readDNAStringSet(ARATHA.cds.file)
## get gene position
ARATHA.cds.genepos <- cds2genepos(ARATHA.cds, "ENSEMBL")
}
}
\seealso{
\code{\link[CRBHits]{isoform2longest}}
}
\author{
Kristian K Ullrich
}
