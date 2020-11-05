#' @title make_KaKs_Calculator2
#' @name make_KaKs_Calculator2
#' @description This function tries to build the prerequisite KaKs_calculator2.0 from source code forked within CRBHits
#' @references Wang D, Zhang Y et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. \emph{Genomics Proteomics Bioinformatics.} \bold{8(1)}, 77-80.
#' @export make_KaKs_Calculator2
#' @author Kristian K Ullrich

make_KaKs_Calculator2 <- function(){
  curwd <- getwd()
  kakscalcpath = paste0(find.package("CRBHits"),
                        "/extdata/KaKs_Calculator2.0/src/")
  if(!dir.exists(kakscalcpath)){
    setwd(paste0(find.package("CRBHits"), "/extdata/"))
    system(paste0("tar -xvf KaKs_Calculator2.0.tar.gz"))
    setwd(paste0(find.package("CRBHits"), "/extdata/KaKs_Calculator2.0/src/"))
    system("make clean")
    system("make")
  }
  setwd(curwd)
}
