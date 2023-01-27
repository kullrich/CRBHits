#' @title cdsfile2aafile
#' @name cdsfile2aafile
#' @description This function translates a cds fasta file into an aa fasta file.
#' @param infile cds fasta file [mandatory]
#' @param outfile aa fasta file [mandatory]
#' @param shorten shorten all sequences to multiple of three [default: FALSE]
#' @param frame  indicates the first base of a the first codon [default: 1]
#' @param framelist  supply vector of frames for each entry [default: NULL]
#' @param genetic.code The genetic code to use for the translation of codons
#' into Amino Acid letters [default: NULL]
#' @return aa fasta file
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom MSA2dist cds2aa
#' @seealso \code{\link[CRBHits]{cds2aa}}
#' @examples
#' ## define file path
#' cdsfile <- system.file("fasta", "ath.cds.fasta.gz", package="CRBHits")
#' ## create empty temp file
#' aafile <- tempfile()
#' ## convert input CDS fasta file into AA fasta file
#' cdsfile2aafile(cdsfile, aafile)
#' aa <- Biostrings::readAAStringSet(aafile)
#' aa
#' @export cdsfile2aafile
#' @author Kristian K Ullrich

cdsfile2aafile <- function(infile,
    outfile,
    shorten=FALSE,
    frame=1,
    framelist=NULL,
    genetic.code=NULL
    ){
    cds <- Biostrings::readDNAStringSet(infile)
    Biostrings::writeXStringSet(MSA2dist::cds2aa(
        cds, shorten=shorten, frame=frame,
        framelist=framelist, genetic.code=genetic.code),
        file=outfile)
}
