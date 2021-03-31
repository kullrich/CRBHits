#' @title dnastring2kaks
#' @name dnastring2kaks
#' @description This function calculates Ka/Ks (pN/pS; according to \emph{Li (1993)} or \emph{Nei and Gojobori (1986)} for all combinations of a \code{DNAStringSet}.
#' @param cds \code{DNAStringSet} A [mandatory]
#' @param align specify if sequence combinations should be aligned pairwise by \code{cds2codonaln}. Default assumes sequences are aligned [default: FALSE]
#' @param model specify codon model either "Li" or "NG86" [default: Li]
#' @param threads number of parallel threads [default: 1]
#' @param ... other codon alignment parameters (see \code{\link[CRBHits]{cds2codonaln}})
#' @return A data.frame of \code{KaKs} values
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq pairwiseAlignment
#' @importFrom seqinr kaks
#' @importFrom doMC registerDoMC
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom tidyr %>% as_tibble pivot_longer
#' @importFrom dplyr slice left_join
#' @seealso \code{\link[seqinr]{kaks}} \code{\link[CRBHits]{cds2codonaln}} \code{\link[Biostrings]{pairwiseAlignment}}
#' @references Nei and Gojobori. (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. \emph{Mol. Biol. Evol.}, \bold{3(5)}, 418-426.
#' @references Ganeshan et al. (1997) Human immunodeficiency virus type 1 genetic evolution in children with different rates of development of disease. \emph{J. Virology.} \bold{71(1)}, 663-677.
#' @references Yang et al. (2000) Codon-substitution models for heterogeneous selection pressure at amino acid sites. \emph{Genetics.} \bold{155(1)}, 431-449.
#' @examples
#' ## load example sequence data
#' data("hiv", package="CRBHits")
#' dnastring2kaks(hiv, model = "Li")
#' dnastring2kaks(hiv, model = "NG86")
#' \dontrun{
#' dnastring2kaks(hiv, model = "NG86", threads = 2)
#' }
#' \dontrun{
#' dnastring2kaks(hiv, model = "NG86", align = TRUE, threads = 2)
#' dnastring2kaks(hiv, model = "Li", align = TRUE, substitutionMatrix = "BLOSUM45")
#' }
#' @export dnastring2kaks
#' @author Kristian K Ullrich

dnastring2kaks <- function(cds,
                           align = FALSE,
                           model = "Li",
                           threads = 1,
                           ...){
  if(!model %in% c("Li", "NG86")){stop("Error: either choose model 'Li' or 'NG86'")}
  if(model == "Li" & align == FALSE){
    OUT <- seqinr::kaks(dnastring2aln(cds))
    OUT.ka <- as.matrix(OUT$ka)
    colnames(OUT.ka) <- rownames(OUT.ka) <- gsub(" ", "", colnames(OUT.ka))
    OUT.ks <- as.matrix(OUT$ks)
    colnames(OUT.ks) <- rownames(OUT.ks) <- gsub(" ", "", colnames(OUT.ks))
    OUT.vka <- as.matrix(OUT$vka)
    colnames(OUT.vka) <- rownames(OUT.vka) <- gsub(" ", "", colnames(OUT.vka))
    OUT.vks <- as.matrix(OUT$vks)
    colnames(OUT.vks) <- rownames(OUT.vks) <- gsub(" ", "", colnames(OUT.vks))
    OUT.ka <- OUT.ka %>% tidyr::as_tibble(rownames = "seq1") %>% tidyr::pivot_longer(colnames(OUT.ka), names_to = "seq2", values_to = "ka") %>% dplyr::slice(uptriidx(ncol(OUT.ka)))
    OUT.ks <- OUT.ks %>% tidyr::as_tibble(rownames = "seq1") %>% tidyr::pivot_longer(colnames(OUT.ks), names_to = "seq2", values_to = "ks") %>% dplyr::slice(uptriidx(ncol(OUT.ks)))
    OUT.vka <- OUT.vka %>% tidyr::as_tibble(rownames = "seq1") %>% tidyr::pivot_longer(colnames(OUT.vka), names_to = "seq2", values_to = "vka") %>% dplyr::slice(uptriidx(ncol(OUT.vka)))
    OUT.vks <- OUT.vks %>% tidyr::as_tibble(rownames = "seq1") %>% tidyr::pivot_longer(colnames(OUT.vks), names_to = "seq2", values_to = "vks") %>% dplyr::slice(uptriidx(ncol(OUT.vks)))
    OUT <- OUT.ka %>% dplyr::left_join(OUT.ks) %>% dplyr::left_join(OUT.vka) %>% dplyr::left_join(OUT.vks)
    OUT <- as.data.frame(OUT)
    attr(OUT, "model") <- "Li"
    attr(OUT, "align") <- "FALSE"
    return(OUT)
  }
  if(model == "Li" & align == TRUE){
    doMC::registerDoMC(threads)
    i <- NULL
    j <- NULL
    OUT <- foreach(i = seq(from = 1, to = length(cds) - 1), .combine=rbind) %dopar% {
      foreach(j = seq(from = i + 1, to = length(cds)), .combine=rbind) %do% {
        c(setNames(i, "Comp1"),
          setNames(j, "Comp2"),
          unlist(seqinr::kaks(dnastring2aln(cds2codonaln(cds[i], cds[j], ...)))))
      }
    }
    OUT <- as.data.frame(OUT)
    attr(OUT, "model") <- "Li"
    attr(OUT, "align") <- "TRUE"
    return(OUT)
  }
  if(model == "NG86" & align == FALSE){
    doMC::registerDoMC(threads)
    i <- NULL
    j <- NULL
    codonmat <- dnastring2codonmat(cds)
    OUT <- foreach(i = seq(from = 1, to = ncol(codonmat) - 1), .combine=rbind) %dopar% {
      foreach(j = seq(from = i + 1, to = ncol(codonmat)), .combine=rbind) %do% {
        c(setNames(i, "Comp1"), setNames(j, "Comp2"), codonmat2pnps(codonmat[, c(i, j)]))
      }
    }
    OUT <- as.data.frame(OUT)
    attr(OUT, "model") <- "NG86"
    attr(OUT, "align") <- "FALSE"
    return(OUT)
  }
  if(model == "NG86" & align == TRUE){
    doMC::registerDoMC(threads)
    i <- NULL
    j <- NULL
    OUT <- foreach(i = seq(from = 1, to = length(cds) - 1), .combine=rbind) %dopar% {
      foreach(j = seq(from = i + 1, to = length(cds)), .combine=rbind) %do% {
        c(setNames(i, "Comp1"),
          setNames(j, "Comp2"),
          codonmat2pnps(dnastring2codonmat(cds2codonaln(cds[i], cds[j], ...))))
      }
    }
    OUT <- as.data.frame(OUT)
    attr(OUT, "model") <- "NG86"
    attr(OUT, "align") <- "TRUE"
    return(OUT)
  }
}
