#' @title cdsfile2aafile
#' @name cdsfile2aafile
#' @description This function translates a cds fasta file into an aa fasta file.
#' @param infile cds fasta file [mandatory]
#' @param outfile aa fasta file [mandatory]
#' @return aa fasta file
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @seealso \code{\link[CRBHits]{cds2aa}}
#' @examples
#' ## define file path
#' cdsfile <- system.file("fasta", "ath.cds.fasta.gz", package = "CRBHits")
#' ## create empty temp file
#' aafile <- tempfile()
#' ## convert input CDS fasta file into AA fasta file
#' cdsfile2aafile(cdsfile, aafile)
#' aa <- Biostrings::readAAStringSet(aafile)
#' aa
#' @export cdsfile2aafile
#' @author Kristian K Ullrich

cdsfile2aafile <- function(infile, outfile){
  cds <- Biostrings::readDNAStringSet(infile)
  Biostrings::writeXStringSet(cds2aa(cds), file = outfile)
}
