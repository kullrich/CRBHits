#' @title cdsfile2aafile
#' @name cdsfile2aafile
#' @description This function translates a cds fasta file into an aa fasta file.
#' @param infile cds fasta file [mandatory]
#' @param outfile aa fasta file [mandatory]
#' @return aa fasta file
#' @importFrom Biostrings readDNAStringSet readAAStringSet writeXStringSet
#' @seealso \code{\link[CRBHits]{cds2aa}}
#' @examples
#' cdsfile <- system.file("fasta/ath.cds.fasta.gz", package = "CRBHits")
#' aafile <- tempfile()
#' cdsfile2aafile(cdsfile, aafile)
#' aa <- Biostrings::readAAStringSet(aafile)
#' aa
#' @export cdsfile2aafile
#' @author Kristian K Ullrich

cdsfile2aafile <- function(infile, outfile){
  cds <- Biostrings::readDNAStringSet(infile)
  Biostrings::writeXStringSet(cds2aa(cds), file = outfile)
}
