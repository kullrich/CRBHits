#' @title rbh2kaks
#' @name rbh2kaks
#' @description This function calculates Ka/Ks (dN/dS; accoring to \emph{Li (1993)} or \emph{Yang and Nielson (2000)} for each (conditional-)reciprocal best hit (CRBHit) pair. The names of the \code{rbh} columns must match the names of the corresponding \code{cds1} and \code{cds2} \code{DNAStringSet} vectors.
#' @param rbhpairs (conditional-)reciprocal best hit (CRBHit) pair result (see \code{\link[CRBHits]{cds2rbh}}) [mandatory]
#' @param cds1 cds1 sequences as \code{DNAStringSet} or \code{url} for first crbh pairs column [mandatory]
#' @param cds2 cds2 sequences as \code{DNAStringSet} or \code{url} for second crbh pairs column [mandatory]
#' @param model specify codon model either "Li" or "YN" [default: Li]
#' @param plotDotPlot specify if dotplot should be plotted (mandatory to define \code{gene.position.cds1} and \code{gene.position.cds1}) [default: FALSE]
#' @param gene.position.cds1 specify gene position for cds1 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param gene.position.cds2 specify gene position for cds2 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param tandem.dups.cds1 specify tandem duplicates for cds1 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param tandem.dups.cds2 specify tandem duplicates for cds2 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param threads number of parallel threads [default: 1]
#' @param kakscalcpath specify the PATH to the KaKs_Calculator binaries [default: /extdata/KaKs_Calculator2.0/src/]
#' @param ... other codon alignment parameters (see \code{\link[CRBHits]{cds2codonaln}})
#' @importFrom doMC registerDoMC
#' @importFrom foreach foreach %dopar%
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @seealso \code{\link[CRBHits]{cds2kaks}},
#' \code{\link[CRBHits]{isoform2longest}},
#' \code{\link[CRBHits]{cds2genepos}}
#' @references Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. \emph{J. Mol. Evol.}, \bold{36}, 96-99.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Yang Z and Nielson R. (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. \emph{Mol. Biol. Evol.}, \bold{17(1)}, 32-43.
#' @examples
#' ##load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ##load example CRBHit pairs
#' data("ath_aly_crbh", package="CRBHits")
#' ##only analyse subset of CRBHit pairs
#' ath_aly_crbh$crbh.pairs <- head(ath_aly_crbh$crbh.pairs)
#' ath_aly_crbh.kaks <- rbh2kaks(rbhpairs = ath_aly_crbh,
#'                               cds1 = ath, cds2 = aly, model = "Li")
#' head(ath_aly_crbh.kaks)
#' @export rbh2kaks
#' @author Kristian K Ullrich

rbh2kaks <- function(rbhpairs, cds1, cds2, model = "Li",
                     plotDotPlot = FALSE,
                     gene.position.cds1 = NULL,
                     gene.position.cds2 = NULL,
                     tandem.dups.cds1 = NULL,
                     tandem.dups.cds2 = NULL,
                     threads = 1,
                     kakscalcpath = paste0(find.package("CRBHits"),
                                           "/extdata/KaKs_Calculator2.0/src/"),
                     ...){
  if(attributes(rbhpairs)$CRBHits.class != "crbh"){
    stop("Please obtain rbhpairs via the cds2rbh() or the cdsfile2rbh() function")
  }
  if(!is.null(gene.position.cds1)){
    if(attributes(gene.position.cds1)$CRBHits.class != "genepos"){
      stop("Please obtain gene position via the cds2genepos() function or add a 'genepos' class attribute")
    }
  }
  if(!is.null(gene.position.cds2)){
    if(attributes(gene.position.cds2)$CRBHits.class != "genepos"){
      stop("Please obtain gene position via the cds2genepos() function or add a 'genepos' class attribute")
    }
  }
  if(!is.null(tandem.dups.cds1)){
    if(attributes(tandem.dups.cds1)$CRBHits.class != "tandemdups"){
      stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
    }
  }
  if(!is.null(tandem.dups.cds2)){
    if(attributes(tandem.dups.cds2)$CRBHits.class != "tandemdups"){
      stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
    }
  }
  if(!model %in% c("Li", "YN")){stop("Error: either choose model 'Li' or 'YN'")}
  if(plotDotPlot){
    if(is.null(gene.position.cds1) & is.null(gene.position.cds2)){
      stop("Error: Please specify gene.position.cds1 and gene.position.cds2")
    }
  }
  #internal function to get cds by name
  get_cds_by_name <- function(x, cds){
    return(cds[names(cds)==x])
  }
  if(class(cds1) == "character"){cds1 <- Biostrings::readDNAStringSet(cds1)}
  if(class(cds2) == "character"){cds2 <- Biostrings::readDNAStringSet(cds2)}
  names(cds1) <- stringr::word(names(cds1), 1)
  names(cds2) <- stringr::word(names(cds2), 1)
  doMC::registerDoMC(threads)
  i <- NULL
  rbhpairs.crbh.pairs <- rbhpairs$crbh.pairs
  rbh.kaks <- foreach::foreach(i = seq(from = 1, to = dim(rbhpairs.crbh.pairs)[1]), .combine = rbind) %dopar% {
    cds2kaks(get_cds_by_name(rbhpairs.crbh.pairs[i,1], cds1), get_cds_by_name(rbhpairs.crbh.pairs[i,2], cds2), model = model, kakscalcpath = kakscalcpath, ...)
  }
  out <- cbind(rbhpairs.crbh.pairs, rbh.kaks)
  attr(out, "CRBHits.class") <- "kaks"
  if(model == "Li"){
    attr(out, "CRBHits.model") <- "Li"
  }
  if(model == "YN"){
    attr(out, "CRBHits.model") <- "YN"
  }
  return(out)
}
