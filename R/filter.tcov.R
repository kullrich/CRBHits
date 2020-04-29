#' @title filter.tcov
#' @name filter.tcov
#' @description This function filters BLAST-like tabular output according to percentage target coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param tcov target coverage [default: 0.0]
#' @return rbh matrix
#' @examples
#' ##load sequence data
#' @export filter.tcov
#' @author Kristian K Ullrich

filter.tcov <- function(rbh, tcov = 0.0){
  return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,14])) >= tcov, , drop = FALSE])
}
