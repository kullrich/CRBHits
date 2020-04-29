#' @title rbh2kaks
#' @name rbh2kaks
#' @description This function calculates ka/ks (dN/dS; accoring to \emph{Li (1993)} or \emph{Yang and Nielson (2000)} for each reciprocal best hit pair. The names of the \code{rbh} columns must match the names of the corresponding \code{cds1} and \code{cds2} \code{DNAStringSet} vectors.
#' @param rbhpairs (conditional-)recirpocal best hit pair matrix [mandatory]
#' @param cds1 cds1 sequences as \code{DNAStringSet} for first crbhpairs column [mandatory]
#' @param cds2 cds2 sequences as \code{DNAStringSet} for second crbhpairs column [mandatory]
#' @param model specify codon model either "Li" or "YN" [default: Li]
#' @param kakscalcpath specify the PATH to the KaKs_Calculator binaries [default: /extdata/KaKs_Calculator2.0/src/]
#' @param threads number of parallel threads [default: 1]
#' @import foreach
#' @import doMC
#' @importFrom seqinr kaks
#' @seealso \code{\link[CRBHits]{cds2kaks}}
#' @references Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. \emph{J. Mol. Evol.}, \bold{36}, 96-99.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Yang Z and Nielson R. (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. \emph{Mol. Biol. Evol.}, \bold{17(1)}, 32-43.
#' @examples
#' ##load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ##load example conditional-reciprocal best hit pair results
#' data("ath_aly_crbh", package="CRBHits")
#' rbh.pairs <- ath_aly_crbh$crbh.pairs
#' ath_aly_crbh.kaks <- rbh2kaks(rbh.pairs, ath, aly)
#' head(ath_aly_crbh.kaks)
#' @export rbh2kaks
#' @author Kristian K Ullrich

rbh2kaks <- function(rbhpairs, cds1, cds2, model = "Li", threads = 1,
                     kakscalcpath = paste0(find.package("CRBHits"),
                                           "/extdata/KaKs_Calculator2.0/src/")){
  if(!model %in% c("Li", "YN")){stop("Error: either choose model 'Li' or 'YN'")}
  #internal function to get cds by name
  get_cds_by_name <- function(x, cds){
    return(cds[names(cds)==x])
  }
  names(cds1) <- unlist(lapply(strsplit(names(cds1), " "), function(x) x[1]))
  names(cds2) <- unlist(lapply(strsplit(names(cds2), " "), function(x) x[1]))
  doMC::registerDoMC(threads)
  i <- NULL
  rbh.kaks <- foreach::foreach(i = seq(from=1, to=dim(rbhpairs)[1]), .combine = rbind) %dopar% {
    cds2kaks(get_cds_by_name(rbhpairs[i,1], cds1), get_cds_by_name(rbhpairs[i,2], cds2), model = model, kakscalcpath = kakscalcpath)
  }
  return(cbind(rbhpairs, rbh.kaks))
}
