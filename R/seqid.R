#' @title seqid
#' @name seqid
#' @description This function returns the first whitespace separated names value of a \code{XStringSet}.
#' @param x \code{XStringSet} [mandatory]
#' @return \code{chr}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @examples
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' seqid(ath)
#' @export seqid
#' @author Kristian K Ullrich

seqid <- function(x){
  if(is.null(names(x))){stop("Error: input needs to have names")}
  if(!is.null(names(x))){
    out <- stringr::word(names(x), 1)
  }
  return(out)
}
