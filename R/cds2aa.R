#' @title cds2aa
#' @name cds2aa
#' @description This function translates a \code{DNAStringSet} into an \code{AAStringSet}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @return \code{AAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq translate
#' @importFrom stringr word
#' @seealso \code{\link[Biostrings]{XStringSet-class}},
#' \code{\link[seqinr]{translate}}
#' @examples
#' ##load example sequence data
#' data("ath", package="CRBHits")
#' cds2aa(ath)
#' @export cds2aa
#' @author Kristian K Ullrich

cds2aa <- function(cds){
  if(class(cds)!="DNAStringSet"){stop("Error: input needs to be a DNAStringSet")}
  if(!is.null(names(cds))){
    names(cds) <- stringr::word(names(cds), 1)
  }
  cds_not_multiple_of_three.idx <- which(Biostrings::width(cds) %% 3 != 0)
  if(length(cds_not_multiple_of_three.idx) > 0){
    cds_not_multiple_of_three <- cds[cds_not_multiple_of_three.idx]
    cds <- cds[-cds_not_multiple_of_three.idx]
  }
  #aa <- Biostrings::AAStringSet(unlist(lapply(as.character(cds), function(x) paste0(seqinr::translate(unlist(strsplit(x, ""))), collapse=""))))
  #aa <- Biostrings::AAStringSet(unlist(lapply(as.character(cds), function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x))))))
  aa <- Biostrings::translate(cds, if.fuzzy.codon = "X")
  return(aa)
}
