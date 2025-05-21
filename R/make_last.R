#' @title make_last
#' @name make_last
#' @description This function tries to build the prerequisite last-1639 from
#' source code forked within CRBHits
#' @param build_dir directory to build last [mandatory]
#' @return compile last and return binary path
#' @importFrom utils unzip
#' @references Kiełbasa SM et al. (2011) Adaptive seeds tame genomic sequence
#' comparison. \bold{Genome Res.} \bold{21} \bold{(3)}, 487-93.
#' @export make_last
#' @author Kristian K Ullrich

make_last <- function(
    build_dir=file.path(find.package("CRBHits"), "extdata")){
    arch <- R.version[["arch"]]
    sysname <- Sys.info()[["sysname"]]
    last_dir <- file.path(build_dir, "last-1639")
    last_path <- file.path(build_dir, "last-1639", "bin")
    lastdb_path <- file.path(build_dir, "last-1639", "bin", "lastdb")
    last_zip <- system.file("extdata",
        "last-1639.zip", package="CRBHits")
    if (!dir.exists(last_dir)) {
        utils::unzip(last_zip, exdir=build_dir)
        Sys.chmod(file.path(last_dir, "build", "gc-inc.sh"), mode="0755")
        Sys.chmod(file.path(last_dir, "build", "mat-doc.sh"), mode="0755")
        Sys.chmod(file.path(last_dir, "build", "mat-inc.sh"), mode="0755")
        Sys.chmod(file.path(last_dir, "build", "seed-doc.sh"), mode="0755")
        Sys.chmod(file.path(last_dir, "build", "seed-inc.sh"), mode="0755")
    }
    if(!file.exists(lastdb_path)){
        # see more installation information here, if make fails
        # last-install-help: http://last.cbrc.jp/doc/last.html
        make_path <- Sys.which("make")
        if (make_path=="") {
            stop("GNU make not found.
                Please ensure Rtools is installed and on your PATH.")
        }
        if (arch %in% c("arm64", "aarch64")) {
            message("Detected ARM architecture (", arch,
                "), patching makefile...")
            makefile1_path <- file.path(last_dir, "makefile")
            if (file.exists(makefile1_path)) {
                makefile1_lines <- readLines(makefile1_path)
                pattern <- "CXXFLAGS = -msse4 -O3 -std=c++11 -pthread"
                repl <- "CXXFLAGS = -O3 -std=c++11 -pthread"
                if (any(grepl(pattern, makefile1_lines, fixed=TRUE))) {
                    message("Patching makefile at: ", makefile1_path)
                    file.copy(makefile1_path, paste0(makefile1_path, ".bak"),
                        overwrite=TRUE)
                    makefile1_lines <- gsub(pattern, repl,
                        makefile1_lines, fixed=TRUE)
                    writeLines(makefile1_lines, makefile1_path)
                } else {
                    message("Expected CXXFLAGS line not found",
                        "- skipping patch.")
                }
            } else {
                warning("Makefile not found at: ", makefile1_path)
            }
            makefile2_path <- file.path(last_dir, "src", "makefile")
            if (file.exists(makefile2_path)) {
              makefile2_lines <- readLines(makefile2_path)
              pattern <- "CXXFLAGS = -msse4 -O3 -Wall -g -std=c++11 -pthread"
              repl <- "CXXFLAGS = -O3 -Wall -g -std=c++11 -pthread"
              if (any(grepl(pattern, makefile2_lines, fixed=TRUE))) {
                message("Patching makefile at: ", makefile2_path)
                file.copy(makefile2_path, paste0(makefile2_path, ".bak"),
                          overwrite=TRUE)
                makefile2_lines <- gsub(pattern, repl,
                    makefile2_lines, fixed=TRUE)
                writeLines(makefile2_lines, makefile2_path)
              } else {
                message("Expected CXXFLAGS line not found",
                        "- skipping patch.")
              }
            } else {
              warning("Makefile not found at: ", makefile2_path)
            }
        } else {
            message("Architecture is ", arch, " - no makefile patch needed.")
        }
        curwd <- getwd()
        on.exit(setwd(curwd))
        setwd(last_dir)
        system2("make")
    }
    return(last_path)
}
