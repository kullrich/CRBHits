#' @title make_dagchainer
#' @name make_dagchainer
#' @description This function tries to build the DAGchainer from source code
#' forked within CRBHits
#' @param build_dir directory to build DAGchainer [mandatory]
#' @return compile DAGchainer and return binary path
#' @importFrom utils unzip
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental
#' genome duplications and synteny. \bold{Bioinformatics}
#' \bold{20} \bold{(18)}, 3643-6.
#' @export make_dagchainer
#' @author Kristian K Ullrich

make_dagchainer <- function(
    build_dir=file.path(find.package("CRBHits"), "extdata")){
    arch <- R.version[["arch"]]
    sysname <- Sys.info()[["sysname"]]
    dagchainer_dir <- file.path(build_dir, "dagchainer")
    dagchainer_path <- file.path(build_dir, "dagchainer", "dagchainer")
    dagchainer_zip <- system.file("extdata",
        "dagchainer.zip", package="CRBHits")
    if (!dir.exists(dagchainer_dir)) {
        utils::unzip(dagchainer_zip, exdir=build_dir)
        Sys.chmod(file.path(dagchainer_dir, "run_DAG_chainer.pl"), mode="0755")
    }
    if(!file.exists(dagchainer_path)){
        make_path <- Sys.which("make")
        if (make_path=="") {
            stop("GNU make not found.
            Please ensure Rtools is installed and on your PATH.")
        }
        curwd <- getwd()
        on.exit(setwd(curwd))
        setwd(dagchainer_dir)
        system2("make")
    }
    return(dagchainer_path)
}
