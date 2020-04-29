#' @title filter.tcov
#' @name filter.tcov
#' @description This function filters BLAST-like tabular output according to percentage target coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param tcov target coverage [default: 0.0]
#' @return rbh matrix
#' @examples
#' ##load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter.tcov(ath_aly_crbh$crbh1, tcov = 0.75))
#' @export filter.tcov
#' @author Kristian K Ullrich

filter.tcov <- function(rbh, tcov = 0.0){
  return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,14])) >= tcov, , drop = FALSE])
}
