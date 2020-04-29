#' @title dnastring2aln
#' @name dnastring2aln
#' @description This function converts a \code{DNAStringSet} into an \code{seqinr} \code{alignment}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @return An object of class \code{alignemnt} which is a list with the following components:
#' \code{nb} the number of aligned sequences
#' \code{nam} a vector of strings containing the names of the aligned sequences
#' \code{seq} a vector of strings containing the aligned sequences
#' \code{com} a vector of strings containing the commentaries for each sequence or \code{NA} if there are no comments
#' @importFrom Biostrings DNAStringSet
#' @importFrom seqinr as.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' @export dnastring2aln
#' @author Kristian K Ullrich

dnastring2aln <- function(dna){
  alignment.nb <- length(dna)
  alignment.nam <- names(dna)
  alignment.seq <- tolower(as.character(dna))
  names(alignment.seq) <- NULL
  alignment.com <- NA
  alignment <- list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
  names(alignment) <- c("nb", "nam", "seq", "com")
  attr(alignment, "class") <- "alignment"
  return(alignment)
}
