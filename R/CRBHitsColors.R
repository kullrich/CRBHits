#' @title CRBHitsColors
#' @name CRBHitsColors
#' @description CRBHits specific color palette.
#' @param n The number of colors (≥ 1) to be in the palette [mandatory]
#' @param alpha.perc an alpha-transparency level in percent [default: 0]
#' @return CRBHits specific color palette
#' @importFrom grDevices col2rgb palette rgb colorRampPalette
#' @seealso \code{\link[grDevices]{palette}}
#' @examples
#' plot(1:9, pch = 20, col = CRBHitsColors(9), cex = 3)
#' ## use alpha
#' plot(1:9, pch = 20, col = CRBHitsColors(9, 50), cex = 3)
#' @export CRBHitsColors
#' @author Kristian K Ullrich

CRBHitsColors <- function(n, alpha.perc = 0){
  my.palette <- palette(c("#CBC106", "#27993C", "#1C6838",
                           "#8EBCB5", "#389CA7", "#4D83AB",
                           "#CB7B26", "#BF565D", "#9E163C"))
  my.n.palette <- colorRampPalette(my.palette)(n)
  my.out.palette <- col2transparent(my.n.palette, alpha.perc)
  return(my.out.palette)
}
