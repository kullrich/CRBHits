#' @title globalDeletion
#' @name globalDeletion
#' @description This function returns a \code{DNAStringSet} reduced by all
#' sites containing any gaps ("-", "+", ".") or missing ("N") sites.
#' @importFrom Biostrings consensusMatrix
#' @param dna \code{DNAStringSet}
#' @export globalDeletion
#' @author Kristian K Ullrich
globalDeletion<-function(dna){
  cM <- Biostrings::consensusMatrix(dna)
  globalDeletionSites <- which(apply(cM, 2, function(x) sum(x[15:18]) >= 1))
  if(length(globalDeletionSites) == 0){
    return(dna)
  }
  return(dnabin2dnastring(dnastring2dnabin(new)[, -globalDeletionSites]))
}
