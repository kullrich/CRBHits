#' @title make.dagchainer
#' @name make.dagchainer
#' @description This function tries to build the DAGchainer from source code forked within CRBHits
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental genome duplications and synteny. \bold{Bioinformatics} \bold{20} \bold{(18)}, 3643-6.
#' @export make.dagchainer
#' @author Kristian K Ullrich

make.dagchainer <- function(){
  curwd <- getwd()
  dagchainerpath <- paste0(find.package("CRBHits"),
         "/extdata/dagchainer/")
  if(!dir.exists(dagchainerpath)){
    setwd(paste0(find.package("CRBHits"), "/extdata/"))
    system(paste0("unzip dagchainer.zip"))
    setwd(paste0(find.package("CRBHits"), "/extdata/dagchainer/"))
  }
  if(!file.exists(paste0(dagchainerpath, "dagchainer"))){
    # see more installation information here, if make fails
    # DAGchainer: http://dagchainer.sourceforge.net/
    system("make")
  }
  setwd(curwd)
}
