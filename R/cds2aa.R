#' @title cds2aa
#' @name cds2aa
#' @description This function translates a \code{DNAStringSet} into an \code{AAStringSet}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @return \code{AAStringSet}
#' @importFrom Biostrings DNAStringSet AAStringSet width
#' @importFrom seqinr translate
#' @seealso \code{\link[Biostrings]{DNAStringSet}},
#' \code{\link[Biostrings]{AAStringSet}},
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
    names(cds) <- unlist(lapply(strsplit(names(cds), " "), function(x) x[1]))
  }
  cds_not_multiple_of_three.idx <- which(Biostrings::width(cds) %% 3 != 0)
  if(length(cds_not_multiple_of_three.idx) > 0){
    cds_not_multiple_of_three <- cds[cds_not_multiple_of_three.idx]
    cds <- cds[-cds_not_multiple_of_three.idx]
  }
  aa <- Biostrings::AAStringSet(unlist(lapply(as.character(cds), function(x) paste0(seqinr::translate(unlist(strsplit(x, ""))), collapse=""))))
  return(aa)
}
