#' @title dnastring2codonmat
#' @name dnastring2codonmat
#' @description This function converts a \code{DNAStringSet} into an \code{codon matrix}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @return An object of class \code{alignment} which is a list with the following components:\cr
#' \code{nb} the number of aligned sequences\cr
#' \code{nam} a vector of strings containing the names of the aligned sequences\cr
#' \code{seq} a vector of strings containing the aligned sequences\cr
#' \code{com} a vector of strings containing the commentaries for each sequence or \code{NA} if there are no comments
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @seealso \code{\link[seqinr]{as.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds1.cds2.aln <- cds2codonaln(cds1, cds2)
#' ## convert into alignment
#' dnastring2codonmat(cds1.cds2.aln)
#' @export dnastring2codonmat
#' @author Kristian K Ullrich

dnastring2codonmat <- function(cds){
  if(class(cds)!="DNAStringSet"){stop("Error: input needs to be a DNAStringSet")}
  if(!is.null(names(cds))){
    names(cds) <- stringr::word(names(cds), 1)
  }
  cds_not_multiple_of_three.idx <- which(Biostrings::width(cds) %% 3 != 0)
  if(length(cds_not_multiple_of_three.idx) > 0){
    cds_not_multiple_of_three <- cds[cds_not_multiple_of_three.idx]
    cds <- cds[-cds_not_multiple_of_three.idx]
  }
  cds.codonmat <- apply(cbind(Biostrings::width(cds), as.character(cds)), 1,
                        function(x) {stringi::stri_sub(x[2], seq(1, as.numeric(x[1]), by = 3),
                                                       length = 3)})
  return(cds.codonmat)
}