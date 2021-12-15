#' @title make_dagchainer
#' @name make_dagchainer
#' @description This function tries to build the DAGchainer from source code
#' forked within CRBHits
#' @return compile DAGchainer
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental
#' genome duplications and synteny. \bold{Bioinformatics}
#' \bold{20} \bold{(18)}, 3643-6.
#' @export make_dagchainer
#' @author Kristian K Ullrich

make_dagchainer <- function(){
    curwd <- getwd()
    dagchainerpath <- paste0(find.package("CRBHits"),
        "/extdata/dagchainer/")
    if(!dir.exists(dagchainerpath)){
        setwd(paste0(find.package("CRBHits"), "/extdata/"))
        system2(command="unzip", "dagchainer.zip")
        setwd(paste0(find.package("CRBHits"), "/extdata/dagchainer/"))
    }
    if(!file.exists(paste0(dagchainerpath, "dagchainer"))){
        # see more installation information here, if make fails
        # DAGchainer: http://dagchainer.sourceforge.net/
        system2("make")
    }
    setwd(curwd)
}
