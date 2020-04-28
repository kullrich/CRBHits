#' @title rbh2kaks.YN
#' @name rbh2kaks.YN
#' @description This function calculates Ka and Ks (dN and dS; accoring to \code{Yang and Nielson 2000}) for each reciprocal best hit pair. The names of the \code{rbh} columns must match the names of the corresponding \code{cds1} and \code{cds2} \code{DNAStringSet} vectors.
#' @param rbh conditional-recirpocal best hit pair matrix [mandatory]
#' @param cds1 cds1 sequences as \code{DNAStringSet} for first rbh column [mandatory]
#' @param cds2 cds2 sequences as \code{DNAStringSet} for second rbh column [mandatory]
#' @param kakscalcpath specify the PATH to the KaKs_Calculator binaries [default: /extdata/KaKs_Calculator2.0/src/]
#' @param threads number of parallel threads [default: 1]
#' @import foreach
#' @import doMC
#' @importFrom seqinr kaks
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Yang Z and Nielson R. (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. \emph{Mol. Biol. Evol.}, \bold{17(1)}, 32-43.
#' @examples
#' ##load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ##load example conditional-reciprocal best hit pair results
#' data("ath_aly_crbh", package="CRBHits")
#' ath_aly_crbh.kaks.YN <- rbh2kaks.YN(ath_aly_crbh$crbh.pairs, ath, aly)
#' head(ath_aly_crbh.kaks.YN)
#' @export rbh2kaks.YN
#' @author Kristian K Ullrich

rbh2kaks.YN <- function(rbh, cds1, cds2,
                        kakscalcpath = paste0(find.package("CRBHits"),
                                              "/extdata/KaKs_Calculator2.0/src/"),
                        threads = 1){
  #internal function to get cds by name
  get_cds_by_name <- function(x, cds){
    return(cds[names(cds)==x])
  }
  if(!dir.exists(kakscalcpath)){stop("Error: KaKs_Calculator2.0 PATH does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.KaKs_Calculator2() function.")}
  if(!file.exists(paste0(kakscalcpath, "AXTConvertor"))){stop("Error: AXTConvertor binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.KaKs_Calculator2() function.")}
  if(!file.exists(paste0(kakscalcpath, "KaKs_Calculator"))){stop("Error: KaKs_Calculator binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.KaKs_Calculator2() function.")}
  names(cds1) <- unlist(lapply(strsplit(names(cds1), " "), function(x) x[1]))
  names(cds2) <- unlist(lapply(strsplit(names(cds2), " "), function(x) x[1]))
  doMC::registerDoMC(threads)
  i <- NULL
  rbh.kaks.YN <- foreach::foreach(i = 1:dim(rbh)[1], .combine = rbind) %dopar% {
    tmp <- tempfile()
    #Biostrings::write.phylip(Biostrings::DNAMultipleAlignment(cds2codonaln(get_cds_by_name(rbh[i,1], cds1), get_cds_by_name(rbh[i,2], cds2))), file = tmp)
    #system(paste0(kakscalcpath, "AXTConvertor ", tmp, " ", tmp, ".axt"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    tmp.codonaln <- cds2codonaln(get_cds_by_name(rbh[i, 1], cds1), get_cds_by_name(rbh[i, 2], cds2))
    sink(tmp)
    cat(paste0(names(tmp.codonaln), collapse="&"), "\n", sep="")
    cat(as.character(tmp.codonaln[1][[1]]), "\n", sep="")
    cat(as.character(tmp.codonaln[2][[1]]), "\n", sep="")
    sink(NULL)
    system(paste0(kakscalcpath, "KaKs_Calculator -i ", tmp, " -o ", tmp, ".YN -m YN"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    out <- unlist(strsplit(readLines(paste0(tmp, ".YN"))[2], "\t"))
    system(paste0("rm ", tmp))
    system(paste0("rm ", tmp, ".YN"))
    #system(paste0("rm ", tmp, ".axt"))
    out
  }
  return(cbind(rbh, rbh.kaks.YN))
}
