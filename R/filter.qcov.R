#' @title filter.qcov
#' @name filter.qcov
#' @description This function filters BLAST-like tabular output according to percentage query coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param qcov query coverage [default: 0.0]
#' @return rbh matrix
#' @examples
#' ##load sequence data
#' @export filter.qcov
#' @author Kristian K Ullrich

filter.qcov <- function(rbh, qcov = 0.0){
  return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,13])) >= qcov, , drop = FALSE])
}
