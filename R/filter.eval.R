#' @title filter.eval
#' @name filter.eval
#' @description This function filters BLAST-like tabular output according to evalue.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param evalue evalue [default: 1e-3]
#' @return rbh matrix
#' @examples
#' ##load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter.eval(ath_aly_crbh$crbh1, evalue = 1e-100))
#' @export filter.eval
#' @author Kristian K Ullrich

filter.eval <- function(rbh, evalue = 1e-3){
  return(rbh[as.numeric(rbh[,11]) <= evalue, , drop = FALSE])
}
