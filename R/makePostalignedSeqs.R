#' @title makePostalignedSeqs
#' @name makePostalignedSeqs
#' @description This function is a fork from an internal function from \code{Biostrings}
#' @param x x
#' @seealso \code{\link[Biostrings]{pairwiseAlignment}},
#' \code{\link[Biostrings]{substitution.matrices}}
#' @export makePostalignedSeqs
#' @author Kristian K Ullrich
makePostalignedSeqs <- get('.makePostalignedSeqs', envir = asNamespace('Biostrings'), inherits = FALSE)
