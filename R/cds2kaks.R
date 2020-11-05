#' @title cds2kaks
#' @name cds2kaks
#' @description This function calculates Ka/Ks (dN/dS; accoring to \emph{Li (1993)} or \emph{Yang and Nielson (2000)} given two coding sequences \code{cds1} and \code{cds2} as \code{DNAStringSet} or \code{DNAString}.
#' @param cds1 single cds1 sequence as \code{DNAStringSet} or \code{DNAString} [mandatory]
#' @param cds2 single cds2 sequence as \code{DNAStringSet} or \code{DNAString} [mandatory]
#' @param model specify codon model either "Li" or "YN" [default: Li]
#' @param kakscalcpath specify the PATH to the KaKs_Calculator binaries [default: /extdata/KaKs_Calculator2.0/src/]
#' @param ... other codon alignment parameters (see \code{\link[CRBHits]{cds2codonaln}})
#' @return vector of ka/ks values as specified by \code{seqinr::kaks} or \code{KaKs_Calculator2.0}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq pairwiseAlignment
#' @importFrom seqinr kaks
#' @seealso \code{\link[seqinr]{kaks}}
#' @references Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. \emph{J. Mol. Evol.}, \bold{36}, 96-99.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Yang Z and Nielson R. (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. \emph{Mol. Biol. Evol.}, \bold{17(1)}, 32-43.
#' @examples
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ## select a sequence pair according to a best hit pair (done for you)
#' cds1 <- ath[1]
#' cds2 <- aly[282]
#' cds2kaks(cds1, cds2, model = "Li")
#' ## using an alternative substitutionMatrix
#' cds2kaks(cds1, cds2, model = "Li", substitutionMatrix = "BLOSUM45")
#' cds2kaks(cds1, cds2, model = "Li", substitutionMatrix = "BLOSUM80")
#' \dontrun{
#' cds2kaks(cds1, cds2, model = "YN")
#' }
#' @export cds2kaks
#' @author Kristian K Ullrich

cds2kaks <- function(cds1, cds2, model = "Li",
                     kakscalcpath = paste0(find.package("CRBHits"),
                                           "/extdata/KaKs_Calculator2.0/src/"),
                     ...){
  if(model == "YN"){
    if(!dir.exists(kakscalcpath)){stop("Error: KaKs_Calculator2.0 PATH does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_KaKs_Calculator2() function.")}
    if(!file.exists(paste0(kakscalcpath, "AXTConvertor"))){stop("Error: AXTConvertor binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_KaKs_Calculator2() function.")}
    if(!file.exists(paste0(kakscalcpath, "KaKs_Calculator"))){stop("Error: KaKs_Calculator binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_KaKs_Calculator2() function.")}
  }
  if(class(cds1) == "DNAStringSet" & length(cds1) != 1){stop("Error: cds1 needs to contain only one sequence")}
  if(class(cds2) == "DNAStringSet" & length(cds2) != 1){stop("Error: cds2 needs to contain only one sequence")}
  if(!model %in% c("Li", "YN")){stop("Error: either choose model 'Li' or 'YN'")}
  if(model == "Li"){
    cds1.cds2.kaks <- unlist(seqinr::kaks(dnastring2aln(cds2codonaln(cds1, cds2, ...))))
    return(cds1.cds2.kaks)
  }
  if(model == "YN"){
    tmp <- tempfile()
    tmp.codonaln <- cds2codonaln(cds1, cds2, ...)
    #create AXT file for KaKs_Calculator2.0
    sink(tmp)
    cat(paste0(names(tmp.codonaln), collapse="&"), "\n", sep="")
    cat(as.character(tmp.codonaln[1][[1]]), "\n", sep="")
    cat(as.character(tmp.codonaln[2][[1]]), "\n", sep="")
    sink(NULL)
    system(paste0(kakscalcpath, "KaKs_Calculator -i ", tmp, " -o ", tmp, ".YN -m YN"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    cds1.cds2.kaks <- unlist(strsplit(readLines(paste0(tmp, ".YN"))[2], "\t"))
    names(cds1.cds2.kaks) <- c("Sequence", "Method", "ka", "ks", "ka/ks", "pValue",
                 "Length", "SSites", "NSites", "FoldSites",  "Substitutions",
                 "SSubstitutions", "NSubstitutions", "FoldSSubstitutions",
                 "FoldNSubstitutions", "DivergenceTime", "SubstitutionRateRatio",
                 "GC", "MLScore", "AICc", "AkaikeWeight", "Model")
    system(paste0("rm ", tmp))
    system(paste0("rm ", tmp, ".YN"))
    return(cds1.cds2.kaks)
  }
}
