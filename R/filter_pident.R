#' @title filter_pident
#' @name filter_pident
#' @description This function filters BLAST-like tabular output according to
#' protein identity.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param pident percent identity [default: 0.0]
#' @param inverse specify if filter should keep the removed values
#' [default: FALSE]
#' @return rbh matrix
#' @importFrom dplyr filter
#' @examples
#' ## load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter_pident(
#'     rbh=ath_aly_crbh$crbh1,
#'     pident=75))
#' @export filter_pident
#' @author Kristian K Ullrich

filter_pident <- function(rbh,
    pident=0.0,
    inverse=FALSE
    ){
    perc_identity <- NULL
    if(inverse){
        #return(rbh[as.numeric(rbh[,3])<pident, , drop=FALSE])
        return(dplyr::filter(rbh, perc_identity<pident))
    } else {
        #return(rbh[as.numeric(rbh[,3])>=pident, , drop=FALSE])
        return(dplyr::filter(rbh, perc_identity>=pident))
    }
}
