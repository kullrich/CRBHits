#' @title filter_alnlen
#' @name filter_alnlen
#' @description This function filters BLAST-like tabular output according to alignment length.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param alnlen alignment length [default: 0.0]
#' @param inverse specify if filter should keep the removed values [default: FALSE]
#' @return rbh matrix
#' @importFrom dplyr filter
#' @examples
#' ## load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter_alnlen(ath_aly_crbh$crbh1, alnlen = 75))
#' @export filter_alnlen
#' @author Kristian K Ullrich

filter_alnlen <- function(rbh, alnlen = 0.0, inverse = FALSE){
  alignment_length <- NULL
  if(inverse){
    #return(rbh[as.numeric(rbh[,4]) < alnlen, , drop = FALSE])
    return(dplyr::filter(rbh, alignment_length < alnlen))
  } else {
    #return(rbh[as.numeric(rbh[,4]) >= alnlen, , drop = FALSE])
    return(dplyr::filter(rbh, alignment_length >= alnlen))
  }
}
