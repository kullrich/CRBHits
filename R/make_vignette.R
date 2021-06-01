#' @title make_vignette
#' @name make_vignette
#' @description This function tries to build the prerequisite last-1243 and KaKs_calculator2.0 from source code forked within CRBHits
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental genome duplications and synteny. \bold{Bioinformatics} \bold{20} \bold{(18)}, 3643-6.
#' @export make_vignette
#' @author Kristian K Ullrich

make_vignette <- function(){
  CRBHits_root <- system.file(package = "CRBHits")
  LastTempDir <- tempdir()
  system(paste0("unzip -o ", CRBHits_root, "/extdata/last-1243.zip -d ", LastTempDir))
  system(paste0("cd ", LastTempDir, "/last-1243/; make"))
  KaKsCalcTempDir <- tempdir()
  system(paste0("tar -C ", KaKsCalcTempDir, " -xvf ", CRBHits_root, "/extdata/KaKs_Calculator2.0.tar.gz"))
  system(paste0("cd ", KaKsCalcTempDir, "/KaKs_Calculator2.0/src/; make clean; make"))
  DAGchainerTempDir <- tempdir()
  system(paste0("unzip -o ", CRBHits_root, "/extdata/dagchainer.zip -d ", DAGchainerTempDir))
  system(paste0("cd ", DAGchainerTempDir, "/dagchainer/; make"))
  return(c(paste0(LastTempDir, "/last-1243/bin/"),
           paste0(KaKsCalcTempDir, "/KaKs_Calculator2.0/src/"),
           paste0(DAGchainerTempDir, "/dagchainer/")))
}
