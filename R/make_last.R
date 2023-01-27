#' @title make_last
#' @name make_last
#' @description This function tries to build the prerequisite last-1418 from
#' source code forked within CRBHits
#' @return compile last
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence
#' comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @export make_last
#' @author Kristian K Ullrich

make_last <- function(){
    curwd <- getwd()
    lastpath <- paste0(find.package("CRBHits"),
        "/extdata/last-1418/")
    if(!dir.exists(lastpath)){
        setwd(paste0(find.package("CRBHits"), "/extdata/"))
        system2(command="unzip", args = "last-1418.zip")
        setwd(paste0(find.package("CRBHits"), "/extdata/last-1418/"))
    }
    if(!file.exists(paste0(lastpath, "bin/lastdb"))){
        # see more installation information here, if make fails
        # last-install-help: http://last.cbrc.jp/doc/last.html
        system2("make")
    }
    setwd(curwd)
}
