#' @title make_vignette
#' @name make_vignette
#' @description This function tries to build the prerequisite last-1639,
#' and DAGchainer from source code forked within CRBHits
#' @return path of prerequisites
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence
#' comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit
#' incorporating gamma-series methods and sliding window strategies.
#' \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental
#' genome duplications and synteny. \bold{Bioinformatics}
#' \bold{20} \bold{(18)}, 3643-6.
#' @export make_vignette
#' @author Kristian K Ullrich

make_vignette <- function(){
    arch <- R.version[["arch"]]
    sysname <- Sys.info()[["sysname"]]
    CRBHits_root <- system.file(package="CRBHits")
    LastTempDir <- tempdir()
    last_path <- make_last(LastTempDir)
    DAGchainerTempDir <- tempdir()
    dagchainer_path <- make_dagchainer(DAGchainerTempDir)
    return(c(
        last_path,
        dagchainer_path))
}
