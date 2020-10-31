#' @title plot.kaks
#' @name plot.kaks
#' @description This function plots Ka/Ks results obtained via `rbh2kaks()` function.
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [mandatory]
#' @param dag specify DAGchainer results as obtained via `rbh2dagchainer()` [default: NULL]
#' @param gene.position.cds1 specify gene position for cds1 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param gene.position.cds2 specify gene position for cds2 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param tandem.dups.cds1 specify tandem duplicates for cds1 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param tandem.dups.cds2 specify tandem duplicates for cds2 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param PlotTitle specify Plot title [default: Ka/Ks results]
#' @param PlotType specify Plot type: "h" histogram, "s" scatterplot or d" dotplot [default: h]
#' @param binw specify binwidth (see \code[\link[ggplot2]{geom_histogram}]) [default: 0.05]
#' @param splitByChr specify if plot should be split by chromosome [default: FALSE]
#' @param colorBy specify if Ka/Ks gene pairs should be colored by "rbh_class", dagchainer", "tandemdups" or "none" [default: none]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @param select.chr filter results for chromosome names [default: NULL]
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap geom_histogram
#' @importFrom gridExtra grid.arrange
#' @export plot.kaks
#' @author Kristian K Ullrich

plot.kaks <- function(kaks,
                      dag = NULL,
                      gene.position.cds1 = NULL,
                      gene.position.cds2 = NULL,
                      tandem.dups.cds1 = NULL,
                      tandem.dups.cds2 = NULL,
                      PlotTitle = "Ka/Ks results",
                      PlotType = "h",
                      binw = 0.05,
                      splitByChr = FALSE,
                      colorBy = "none",
                      ka.max = 5,
                      ks.max = 5,
                      ka.min = 0,
                      ks.min = 0,
                      select.chr = NULL){
  if(attributes(kaks)$CRBHits.class != "kaks"){
      stop("")
  }
  if(!is.null(dag)){
    if(attributes(dag)$CRBHits.class != "dagchainer"){
      stop("")
    }
  }
  kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x > ka.max, NA, x)))
  kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x < ka.min, NA, x)))
  kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x > ks.max, NA, x)))
  kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x < ks.min, NA, x)))
  if(attributes(kaks)$CRBHits.model == "Li"){
    kaks <- kaks %>% dplyr::mutate(kaks = kaks$ka / kaks$ks)
    kaks$kaks[is.infinite(kaks$kaks)] <- NA
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
    dagchainer.pair.ids <- apply(cbind(dag$gene1.seq.id, dag$gene2.seq.id),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    kaks.pair.ids <- apply(cbind(kaks$aa1, kaks$aa2),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    kaks.results.dag <- dag$dagchainer_group[match(kaks.pair.ids, dagchainer.pair.ids)]
    kaks <- kaks %>% dplyr::mutate(dag = kaks.results.dag)
    g <- kaks %>% subset(!is.na(ka)) %>%
      subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
      ggplot2::ggplot()
    g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka)) +
      ggplot2::theme(legend.position="none") +
      ggplot2::geom_abline(aes(slope = 1, intercept = 0)) +
      ggplot2::geom_point(data = subset(g$data, !is.na(dag)),
                          shape = 20, aes(x = ks, y = ka, col = dag)) +
      ggplot2::scale_colour_manual(
        values = CRBHitsColors(length(unique(kaks$dag)))) +
      ggplot2::ggtitle(PlotTitle)
    g.ka <- kaks %>% subset(!is.na(ka)) %>% ggplot2::ggplot() +
      ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = dag)) +
      ggplot2::theme(legend.position="none") +
      ggplot2::scale_fill_manual(
        values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
      ggplot2::ggtitle("Ka")
    g.ks <- kaks %>% subset(!is.na(ks)) %>% ggplot2::ggplot() +
      ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = dag)) +
      ggplot2::theme(legend.position="none") +
      ggplot2::scale_fill_manual(
        values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
      ggplot2::ggtitle("Ks")
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
