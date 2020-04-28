#' @title filter.qcov
#' @name filter.qcov
#' @description This function filters BLAST-like tabular output according to percentage query coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param qcov query coverage [default: 0.0]
#' @examples
#' ##load sequence data
#' @export filter.qcov
#' @author Kristian K Ullrich

filter.qcov <- function(rbh, qcov = 0.0){
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[4])/as.numeric(x[13]) >= qcov, TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
