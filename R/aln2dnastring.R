#' @title aln2dnastring
#' @name aln2dnastring
#' @description This function converts a \code{seqinr} \code{alignment} into an \code{DNAStringSet}.
#' @param aln \code{seqinr} \code{alignment} [mandatory]
#' @return An object of class \code{DNAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @seealso \code{\link[seqinr]{as.alignment}} \code{\link[Biostrings]{DNAStringSet}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds1.cds2.aln <- cds2codonaln(cds1, cds2)
#' ## convert into alignment
#' dnastring2aln(cds1.cds2.aln)
#' ## convert into DNAStringSet
#' aln2dnastring(dnastring2aln(cds1.cds2.aln))
#' @export aln2dnastring
#' @author Kristian K Ullrich

aln2dnastring <- function(aln){
  dna <- Biostrings::DNAStringSet(unlist(aln$seq))
  names(dna) <- aln$nam
  return(dna)
}
