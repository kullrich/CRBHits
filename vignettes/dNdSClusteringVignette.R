## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#plot expected pident by alignment length using eq2 from Rost (1999)
get_pident_by_length <- function(x){
  eq2 <- function(L){
    if(L <= 11){return(100)}
    if(L <= 450){return(480*(L^(-0.32*(1+(exp(-L/1000))))))}
    if(L > 450){return(19.5)}
  }
  return(unlist(lapply(x, eq2)))
}
curve(get_pident_by_length, 11, 500, pch = 20, xlab = "alignment length", ylab = "pident",
      main = "expected protein identity (eq2; Rost B. 19999)")

