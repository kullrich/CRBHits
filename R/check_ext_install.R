#' @title check_ext_install
#' @name check_ext_install
#' @description This function checks if external binary exists
#' @param ext_name external tool name [mandatory]
#' @param ext_dir external tool directory [mandatory]
#' @param binary_name external binary name (per default ext_name) [mandatory]
#' @return path of prerequisites
#' @export check_ext_install
#' @author Kristian K Ullrich

check_ext_install <- function(
    ext_name,
    ext_dir,
    binary_name=ext_name){
    if (!dir.exists(ext_dir)){
        stop(sprintf(
            "Error: directory for '%s' not found at '%s'.\n
            Please provide the correct path and ensure installation
            prerequisites are met.",
            ext_name, ext_dir))
    }
    binary_path <- file.path(ext_dir, binary_name)
    if (!file.exists(binary_path)) {
        stop(sprintf(
            "Error: binary '%s' not found in '%s'.\n
            Please provide the correct path or check the tool installation.",
            binary_name, ext_dir))
    }
    return(normalizePath(binary_path))
}
