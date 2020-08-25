#' @title make.vignette
#' @name make.vignette
#' @description This function tries to build the prerequisite last-1080 and KaKs_calculator2.0 from source code forked within CRBHits
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @export make.vignette
#' @author Kristian K Ullrich

make.vignette <- function(){
  CRBHits_root <- system.file(package = "CRBHits")
  LastTempDir <- tempdir()
  system(paste0("unzip -o ", CRBHits_root, "/extdata/last-1080.zip -d ", LastTempDir))
  system(paste0("cd ", LastTempDir, "/last-1080/; make"))
  KaKsCalcTempDir <- tempdir()
  system(paste0("tar -C ", KaKsCalcTempDir, " -xvf ", CRBHits_root, "/extdata/KaKs_Calculator2.0.tar.gz"))
  system(paste0("cd ", KaKsCalcTempDir, "/KaKs_Calculator2.0/src/; make clean; make"))
  return(c(paste0(LastTempDir, "/last-1080/src/"),
           paste0(KaKsCalcTempDir, "/KaKs_Calculator2.0/src/")))
}
