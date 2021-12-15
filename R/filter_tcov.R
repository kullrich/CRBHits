#' @title filter_tcov
#' @name filter_tcov
#' @description This function filters BLAST-like tabular output according to
#' percentage target coverage.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param tcov target coverage [default: 0.0]
#' @param inverse specify if filter should keep the removed values
#' [default: FALSE]
#' @return rbh matrix
#' @importFrom dplyr filter
#' @examples
#' ## load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter_tcov(ath_aly_crbh$crbh1, tcov=0.75))
#' @export filter_tcov
#' @author Kristian K Ullrich

filter_tcov <- function(rbh, tcov=0.0, inverse=FALSE){
    alignment_length <- NULL
    subject_length <- NULL
    if(inverse){
        #return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,14]))<tcov, ,
        #drop=FALSE])
        return(dplyr::filter(rbh, (alignment_length / subject_length)<tcov))
    } else {
          #return(rbh[(as.numeric(rbh[,4])/as.numeric(rbh[,14]))>=tcov, ,
          #drop=FALSE])
          return(dplyr::filter(rbh, (alignment_length / subject_length)>=tcov))
    }
}
