#' @title rbh2kaks
#' @name rbh2kaks
#' @description This function calculates Ka and Ks (dN and dS; accoring to *Li 1993*) for each reciprocal best hit pair. The names of the \code{rbh} columns must match the names of the corresponding \code{cds1} and \code{cds2} \code{DNAStringSet} vectors.
#' @param rbh conditional-recirpocal best hit pair matrix [mandatory]
#' @param cds1 cds1 sequences as \code{DNAStringSet} for first rbh column [mandatory]
#' @param cds2 cds2 sequences as \code{DNAStringSet} for second rbh column [mandatory]
#' @param threads number of parallel threads [default: 1]
#' @import foreach
#' @import doMC
#' @importFrom seqinr kaks
#' @references Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. \emph{J. Mol. Evol.}, \bold{36}, 96-99.
#' @examples
#' ##load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ##load example conditional-reciprocal best hit pair results
#' data("ath_aly_crbh", package="CRBHits")
#' ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh$crbh.pairs, ath, aly)
#' head(ath_aly_crbh.kaks)
#' @export rbh2kaks
#' @author Kristian K Ullrich

rbh2kaks <- function(rbh, cds1, cds2, threads = 1){
  #internal function to get cds by name
  get_cds_by_name <- function(x, cds){
    return(cds[names(cds)==x])
  }
  names(cds1) <- unlist(lapply(strsplit(names(cds1), " "), function(x) x[1]))
  names(cds2) <- unlist(lapply(strsplit(names(cds2), " "), function(x) x[1]))
  doMC::registerDoMC(threads)
  i <- NULL
  rbh.kaks <- foreach::foreach(i = 1:dim(rbh)[1], .combine = rbind) %dopar% {
    unlist(seqinr::kaks(dnastring2aln(cds2codonaln(
      get_cds_by_name(rbh[i,1], cds1),
      get_cds_by_name(rbh[i,2], cds2)))))
  }
  return(cbind(rbh, rbh.kaks))
}

