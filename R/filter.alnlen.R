#' @title filter.alnlen
#' @name filter.alnlen
#' @description This function filters BLAST-like tabular output according to alignment length.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param alnlen alignment length [default: 0.0]
#' @examples
#' ##load sequence data
#' @export filter.alnlen
#' @author Kristian K Ullrich

filter.alnlen <- function(rbh, alnlen = 0.0){
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[4]) >= alnlen, TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
