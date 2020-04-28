#' @title filter.tcov
#' @name filter.tcov
#' @description This function filters BLAST-like tabular output according to percentage target coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param tcov target coverage [default: 0.0]
#' @examples
#' ##load sequence data
#' @export filter.tcov
#' @author Kristian K Ullrich

filter.tcov <- function(rbh, tcov = 0.0){
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[4])/as.numeric(x[14]) >= tcov, TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
