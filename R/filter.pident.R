#' @title filter.pident
#' @name filter.pident
#' @description This function filters BLAST-like tabular output according to protein identity.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param pident percent identity [default: 0.0]
#' @return rbh matrix
#' @examples
#' ##load sequence data
#' @export filter.pident
#' @author Kristian K Ullrich

filter.pident <- function(rbh, pident = 0.0){
  return(rbh[as.numeric(rbh[,3]) >= pident, , drop = FALSE])
}
