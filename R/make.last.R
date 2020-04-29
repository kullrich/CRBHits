#' @title make.last
#' @name make.last
#' @description This function tries to build the prerequisite last-1060 from source code forked within CRBHits
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @export make.last
#' @author Kristian K Ullrich

make.last <- function(){
  curwd <- getwd()
  lastpath <- paste0(find.package("CRBHits"),
         "/extdata/last-1060/src/")
  if(!dir.exists(lastpath)){
    setwd(paste0(find.package("CRBHits"), "/extdata/"))
    system(paste0("unzip last-1060.zip"))
    setwd(paste0(find.package("CRBHits"), "/extdata/last-1060/"))
  }
  if(!file.exists(paste0(lastpath, "lastdb"))){
    # see more installation iformation here, if make fails
    # last-install-help: http://last.cbrc.jp/doc/last.html
    system("make")
  }
  setwd(curwd)
}
