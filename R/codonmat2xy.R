#' @title codonmat2xy
#' @name codonmat2xy
#' @description This function calculates average behavior of each codon for all pairwise comparisons for indels, syn, and nonsyn mutations according to \emph{Nei and Gojobori (1986)}.
#' @param codonmat \code{codonmat} A [mandatory]
#' @param threads number of parallel threads [default: 1]
#' @return A \code{data.frame} object with the following components:\cr
#' \code{Codon} Codon index\cr
#' \code{n} number of comparison\cr
#' \code{SynSum} Sum of syn\cr
#' \code{NonSynSum} Sum of nonsyn\cr
#' \code{IndelSum} Sum of indels\cr
#' \code{SynMean} average syn per codon\cr
#' \code{NonSynMean} average nonsyn per codon\cr
#' \code{IndelMean} average indels per codon\cr
#' \code{CumSumSynMean} cumulative average syn per codon\cr
#' \code{CumSumNonSynMean} cumulative average nonsyn per codon\cr
#' \code{CumSumIndelMean} cumulative indels per codon\cr
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq pairwiseAlignment
#' @importFrom seqinr kaks
#' @importFrom tidyr %>% unite
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom dplyr group_by filter count left_join summarise mutate
#' @seealso \code{\link[seqinr]{kaks}} \code{\link[CRBHits]{codonmat2pnps}} \code{\link[CRBHits]{dnastring2kaks}}
#' @references Nei and Gojobori. (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. \emph{Mol. Biol. Evol.}, \bold{3(5)}, 418-426.
#' @references Ganeshan et al. (1997) Human immunodeficiency virus type 1 genetic evolution in children with different rates of development of disease. \emph{J. Virology.} \bold{71(1)}, 663-677.
#' @references Yang et al. (2000) Codon-substitution models for heterogeneous selection pressure at amino acid sites. \emph{Genetics.} \bold{155(1)}, 431-449.
#' @examples
#' ## load example sequence data
#' data("hiv", package="CRBHits")
#' hiv.xy <- codonmat2xy(dnastring2codonmat(hiv))
#' @export codonmat2xy
#' @author Kristian K Ullrich

codonmat2xy <- function(codonmat, threads = 1){
  doMC::registerDoMC(threads)
  i <- NULL
  j <- NULL
  k <- NULL
  OUT <- foreach(i = seq(from = 1, to = nrow(codonmat)), .combine=rbind) %dopar% {
    foreach(j = seq(from = 1, to = ncol(codonmat) - 1), .combine=rbind) %do% {
      foreach(k = seq(from = j + 1, to = ncol(codonmat)), .combine=rbind) %do% {
        c(setNames(i, "Codon"),
          setNames(j, "Comp1"),
          setNames(k, "Comp2"),
          setNames(compareCodons(codonmat[i, j], codonmat[i, k]), c("syn", "nonsyn", "indel")))
      }
    }
  }
  OUT <- as.data.frame(OUT)
  OUT.NAs <- OUT %>% dplyr::group_by(Codon) %>% dplyr::filter(!is.na(syn)) %>%
                     dplyr::count(Codon)
  OUT.SynSum <- OUT %>% dplyr::group_by(Codon) %>%
                        dplyr::summarise(SynSum = sum(syn, na.rm = TRUE))
  OUT.NonSynSum <- OUT %>% dplyr::group_by(Codon) %>%
                           dplyr::summarise(NonSynSum = sum(nonsyn, na.rm = TRUE))
  OUT.IndelSum <- OUT %>% dplyr::group_by(Codon) %>%
                          dplyr::summarise(IndelSum = sum(indel, na.rm = TRUE))
  OUT.join <- dplyr::left_join(OUT.NAs, OUT.SynSum) %>%
             dplyr::left_join(OUT.NonSynSum) %>% dplyr::left_join(OUT.IndelSum)
  OUT.xy <- OUT.join %>% dplyr::mutate(SynMean = SynSum/n,
                                       NonSynMean = NonSynSum/n,
                                       IndelMean = IndelSum/n)
  OUT.xy <- OUT.xy %>% dplyr::mutate(CumSumSynMean = cumsum(SynMean),
                                       CumSumNonSynMean = cumsum(NonSynMean),
                                       CumSumIndelMean = cumsum(IndelMean))
  return(OUT.xy)
}
