#' @title filter_rost1999
#' @name filter_rost1999
#' @description This function filters BLAST-like tabular output according to equation 2 of Rost 1999.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @param inverse specify if filter should keep the removed values [default: FALSE]
#' @return rbh matrix
#' @references Rost B. (1999). Twilight zone of protein sequence alignments. \emph{Protein Engineering}, \bold{12(2)}, 85-94.
#' @examples
#' ## load crbh data
#' data(ath_aly_crbh)
#' dim(ath_aly_crbh$crbh1)
#' dim(filter_rost1999(ath_aly_crbh$crbh1))
#' @export filter_rost1999
#' @author Kristian K Ullrich

filter_rost1999 <- function(rbh, inverse = FALSE){
  #internal function to calculate pident by length
  get_pident_by_length <- function(x){
    eq2 <- function(x){
      if(x <= 11){return(100)}
      if(x <= 450){return(480*(x^(-0.32*(1+(exp(-x/1000))))))}
      if(x > 450){return(19.5)}
    }
    return(unlist(lapply(x, eq2)))
  }
  pident_by_length <- get_pident_by_length(as.numeric(rbh[,4]))
  if(inverse){
    return(rbh[as.numeric(rbh[,3]) < pident_by_length, , drop = FALSE])
  } else {
    return(rbh[as.numeric(rbh[,3]) >= pident_by_length, , drop = FALSE])
  }
}
