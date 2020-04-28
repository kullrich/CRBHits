#' @title filter.eval
#' @name filter.eval
#' @description This function filters BLAST-like tabular output according to evalue.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param evalue evalue [default: 1e-3]
#' @examples
#' ##load sequence data
#' data("blastTab", package="CRBHits")
#' dim(blastTab)
#' blastTab.filtered <- blastTab[apply(blastTab, 1, function(x){
#' filter.eval(x, 1e-3)
#' }), ,drop=FALSE]
#' @export filter.eval
#' @author Kristian K Ullrich

filter.eval <- function(rbh, evalue = 1e-3){
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[11]) <= evalue, TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
