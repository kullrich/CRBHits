#' @title filter_eval
#' @name filter_eval
#' @description This function filters BLAST-like tabular output according to
#' evalue.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param eval evalue [default: 1e-3]
#' @param inverse specify if filter should keep the removed values
#' [default: FALSE]
#' @return rbh matrix
#' @importFrom dplyr filter
#' @examples
#' ## load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter_eval(
#'     rbh=ath_aly_crbh$crbh1,
#'     eval=1e-100))
#' @export filter_eval
#' @author Kristian K Ullrich

filter_eval <- function(rbh,
    eval=1e-3,
    inverse=FALSE
    ){
    evalue <- NULL
    if(inverse){
        #return(rbh[as.numeric(rbh[,11])>evalue, , drop=FALSE])
        return(dplyr::filter(rbh, evalue>eval))
    } else {
        #return(rbh[as.numeric(rbh[,11])<=evalue, , drop=FALSE])
        return(dplyr::filter(rbh, evalue<=eval))
    }
}
