#' @title filter.qcov
#' @name filter.qcov
#' @description This function filters BLAST-like tabular output according to percentage query coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param qcov query coverage [default: 0.0]
#' @param inverse specify if filter should keep the removed values [default: FALSE]
#' @return rbh matrix
#' @examples
#' ##load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter.qcov(ath_aly_crbh$crbh1, qcov = 0.75))
#' @export filter.qcov
#' @author Kristian K Ullrich

filter.qcov <- function(rbh, qcov = 0.0, inverse = FALSE){
  if(inverse){
    return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,13])) < qcov, , drop = FALSE])
  } else {
    return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,13])) >= qcov, , drop = FALSE])
  }
}
