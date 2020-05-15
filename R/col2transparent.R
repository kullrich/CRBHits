#' @title col2transparent
#' @name col2transparent
#' @description R color to transparent conversion.
#' @param col R color name or R color pal [mandatory]
#' @param alpha.perc an alpha-transparency level in percent [default: 0]
#' @importFrom grDevices col2rgb palette rgb
#' @seealso \code{\link[grDevices]{palette}}
#' @examples
#' col2transparent("green", 50)
#' ## Demonstrate the colors 1:8 in different palettes using a custom matplot()
#' sinplot <- function(main=NULL) {
#'   x <- outer(
#'     seq(-pi, pi, length.out = 50),
#'     seq(0, pi, length.out = 8),
#'     function(x, y) sin(x - y)
#'   )
#'   matplot(x, type = "l", lwd = 4, lty = 1, col = 1:8, ylab = "", main=main)
#' }
#' palette("R3"); sinplot("R3")
#' palette(col2transparent(palette("R3"), 50)); sinplot("R3 - transparent")
#' @export col2transparent
#' @author Kristian K Ullrich

col2transparent <- function(col, alpha.perc = 0){
  if(class(col) == "character"){
    if(length(col) == 1){
      alpha = (100 - alpha.perc) * 255 / 100
      R = col2rgb(col)[1]
      G = col2rgb(col)[2]
      B = col2rgb(col)[3]
      return(rgb(R, G, B, alpha, max = 255))
    } else {
      return(unlist(lapply(col,
                           function(x) {col2transparent(x, alpha.perc)}
      )))
    }
  }
  if(class(col) == "palette"){
    return(palette(unlist(lapply(col,
                                 function(x) {col2transparent(x, alpha.perc)}
                                 ))))
  }
}
