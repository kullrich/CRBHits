#' @title plot.kaks
#' @name plot.kaks
#' @description This function plots Ka/Ks results obtained via `rbh2kaks()` function.
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [mandatory]
#' @param dag specify DAGchainer results as obtained via `rbh2dagchainer()` [default: NULL]
#' @param tandems specify tandem duplicates results as obtained via `tandemdups()` [default: NULL]
#' @param PlotTitle specify Plot title [default: Ka/Ks results]
#' @param colorBy specify if Ka/Ks gene pairs should be colored by "rbh_class", dagchainer", "tandemdups" or "none" [default: none]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap geom_histogram
#' @importFrom gridExtra grid.arrange
#' @export plot.kaks
#' @author Kristian K Ullrich

plot.kaks <- function(kaks, dag = NULL, tandems = NULL, PlotTitle = "Ka/Ks results", colorBy = "none", ka.max = 5, ks.max = 5, ka.min = 0, ks.min = 0){
  if(attributes(kaks)$CRBHits.class != "kaks"){
      stop("")
  }
  if(!is.null(dag)){
    if(attributes(dag)$CRBHits.class != "dagchainer"){
      stop("")
    }
  }
  if(colorBy == "none"){
    g <- kaks %>% ggplot2::ggplot()
    g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka)) +
      ggplot2::ggtitle(PlotTitle)
    g.ka <- g + ggplot2::geom_histogram(aes(x = ka)) + ggplot2::ggtitle("Ka")
    g.ks <- g + ggplot2::geom_histogram(aes(x = ks)) + ggplot2::ggtitle("Ks")
    gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
  }
  if(colorBy == "rbh_class"){
    g <- kaks %>% ggplot2::ggplot()
    g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka)) +
      ggplot2::ggtitle(PlotTitle)
    g.ka <- g + ggplot2::geom_histogram(aes(x = ka)) + ggplot2::ggtitle("Ka")
    g.ks <- g + ggplot2::geom_histogram(aes(x = ks)) + ggplot2::ggtitle("Ks")
    gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
  }
  if(colorBy == "dagchainer"){
    g <- kaks %>% ggplot2::ggplot()
    g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka)) +
      ggplot2::ggtitle(PlotTitle)
    g.ka <- g + ggplot2::geom_histogram(aes(x = ka)) + ggplot2::ggtitle("Ka")
    g.ks <- g + ggplot2::geom_histogram(aes(x = ks)) + ggplot2::ggtitle("Ks")
    gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
  }
  if(colorBy == "tandemdups"){
    g <- kaks %>% ggplot2::ggplot()
    g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka)) +
      ggplot2::ggtitle(PlotTitle)
    g.ka <- g + ggplot2::geom_histogram(aes(x = ka)) + ggplot2::ggtitle("Ka")
    g.ks <- g + ggplot2::geom_histogram(aes(x = ks)) + ggplot2::ggtitle("Ks")
    gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
  }
}
