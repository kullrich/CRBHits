#' @title filter.alnlen
#' @name filter.alnlen
#' @description This function filters BLAST-like tabular output according to alignment length.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param alnlen alignment length [default: 0.0]
#' @return rbh matrix
#' @examples
#' ##load sequence data
#' @export filter.alnlen
#' @author Kristian K Ullrich

filter.alnlen <- function(rbh, alnlen = 0.0){
  return(rbh[as.numeric(rbh[,4]) >= alnlen, , drop = FALSE])
}
