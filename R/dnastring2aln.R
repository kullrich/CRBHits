#' @title dnastring2aln
#' @name dnastring2aln
#' @description This function converts a \code{DNAStringSet} into an \code{seqinr} \code{alignment}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @importFrom Biostrings DNAStringSet
#' @importFrom seqinr as.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' @export dnastring2aln
#' @author Kristian K Ullrich

dnastring2aln <- function(dna){
  alignment.nb <- length(dna)
  alignment.nam <- names(dna)
  alignment.seq <- tolower(as.character(dna))
  alignment.com <- NA
  alignment <- list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
  names(alignment) <- c("nb", "nam", "seq", "com")
  attr(alignment, "class") <- "alignment"
  return(alignment)
}
