#' @title filter.pident
#' @name filter.pident
#' @description This function filters BLAST-like tabular output according to protein identity.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param pident percent identity [default: 0.0]
#' @examples
#' ##load sequence data
#' @export filter.pident
#' @author Kristian K Ullrich

filter.pident <- function(rbh, pident = 0.0){
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[3]) >= pident, TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
