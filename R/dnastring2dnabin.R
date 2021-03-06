#' @title dnastring2dnabin
#' @name dnastring2dnabin
#' @description This function converts a \code{DNAStringSet} into an \code{ape} \code{DNAbin}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @return An object of class \code{DNAbin}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @importFrom ape as.DNAbin.alignment
#' @seealso \code{\link[seqinr]{as.alignment}} \code{\link[ape]{as.DNAbin.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds1.cds2.aln <- cds2codonaln(cds1, cds2)
#' ## convert into alignment
#' dnastring2dnabin(cds1.cds2.aln)
#' @export dnastring2dnabin
#' @author Kristian K Ullrich

dnastring2dnabin <- function(dna){
  if(class(dna)!="DNAStringSet"){stop("Error: input needs to be a DNAStringSet")}
  alignment.nb <- length(dna)
  alignment.nam <- names(dna)
  alignment.seq <- tolower(as.character(dna))
  names(alignment.seq) <- NULL
  alignment.com <- NA
  alignment <- list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
  names(alignment) <- c("nb", "nam", "seq", "com")
  attr(alignment, "class") <- "alignment"
  alignment.bin <- ape::as.DNAbin.alignment(alignment)
  return(alignment.bin)
}
