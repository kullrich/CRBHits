#' @title filter.rost1999
#' @name filter.rost1999
#' @description This function filters BLAST-like tabular output according to equation 2 of Rost 1999.
#' @param rbh BLAST-like tabular matrix [mandatory]
#' @references Rost B. (1999). Twilight zone of protein sequence alignments. \emph{Protein Engineering}, \bold{12(2)}, 85-94.
#' @examples
#' ##load sequence data
#' @export filter.rost1999
#' @author Kristian K Ullrich

filter.rost1999 <- function(rbh){
  #internal function to calculate pident by length
  get_pident_by_length <- function(x){
    if(x<=11){return(100)}
    if(x<=450){return(480*(x^(-0.32*(1+(exp(-x/1000))))))}
    if(x>450){return(19.5)}
  }
  retain.idx <- apply(rbh, 1, function(x) {ifelse(as.numeric(x[3]) >= get_pident_by_length(as.numeric(x[4])), TRUE, FALSE)})
  return(rbh[retain.idx, , drop = FALSE])
}
