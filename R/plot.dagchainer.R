#' @title plot.dagchainer
#' @name plot.dagchainer
#' @description This function plots DAGchainer (http://dagchainer.sourceforge.net/) results obtained via `rbh2dagchainer()` function.
#' @param dag specify DAGchainer results as obtained via `rbh2dagchainer()` [mandatory]
#' @param DotPlotTitle specify DotPlot title [default: DAGchainer results]
#' @param colorBy specify if dagchainer groups should be colored by "Ka", "Ks", "Ka/Ks" or "none" [default: none]
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [default: NULL]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @param select.chr filter results for chromosome names [default: NULL]
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap scale_colour_manual scale_colour_continuous
#' @examples
#' ## load example sequence data
#' data("ath_aly_ncbi_dagchainer", package="CRBHits")
#' ## plot DAGchainer results - default
#' plot.kaks(ath_aly_ncbi_dagchainer)
#' @export plot.dagchainer
#' @author Kristian K Ullrich

plot.dagchainer <- function(dag,
                            DotPlotTitle = "DAGchainer results",
                            colorBy = "none",
                            kaks = NULL,
                            ka.max = 5,
                            ks.max = 5,
                            ka.min = 0,
                            ks.min = 0,
                            select.chr = NULL){
  if(attributes(dag)$CRBHits.class != "dagchainer"){
    stop("Please obtain DAGchainer results via the rbh2dagchainer() function or add a 'dagchainer' class attribute")
  }
  if(!is.null(kaks)){
    if(attributes(kaks)$CRBHits.class != "kaks"){
      stop("Please obtain Ka/Ks via the rbh2kaks() function or add a 'kaks' class attribute")
    }
  }
  selfblast <- attributes(dag)$selfblast
  if(!is.null(select.chr)){
    dag <- dag %>% dplyr::filter(gene1.chr %in% select.chr) %>% dplyr::filter(gene2.chr %in% select.chr)
  }
  if(colorBy == "none"){
    dagchainer.results.split <- dag %>% group_by(gene1.chr, gene2.chr) %>% group_split(.keep = TRUE)
    g <- dag %>% group_by(gene1.chr, gene2.chr) %>% ggplot2::ggplot()
    if(selfblast){
      g.none <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = as.factor(dagchainer_group))) +
        ggplot2::theme(legend.position="none") +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_point(data = subset(g$data, gene1.chr == gene2.chr), shape = 20,
                            aes(x = gene1.mid, y = gene2.mid,
                                col = as.factor(dagchainer_group))) +
        ggplot2::scale_colour_manual(
          values = CRBHitsColors(length(unique(dag$dagchainer_group)))) +
        geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                    aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = CRBHitsColors(length(unique(dag$dagchainer_group)))) +
        ggplot2::ggtitle(DotPlotTitle)
      return(g.none)
    }
    if(!selfblast){
      g.none <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = as.factor(dagchainer_group))) +
        ggplot2::theme(legend.position="none") +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::scale_colour_manual(
          values = CRBHitsColors(length(unique(dag$dagchainer_group)))) +
        ggplot2::ggtitle(DotPlotTitle)
      return(g.none)
    }
  }
  if(colorBy == "Ka" | colorBy == "Ks" | colorBy == "Ka/Ks"){
    if(is.null(kaks)){
      stop("kaks needs to be specified and can be obtained via the rbh2kaks() function")
    }
    dagchainer.pair.ids <- apply(cbind(dag$gene1.seq.id, dag$gene2.seq.id),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    kaks.pair.ids <- apply(cbind(kaks$aa1, kaks$aa2),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    dagchainer.results.ka <- kaks$ka[match(dagchainer.pair.ids, kaks.pair.ids)]
    dagchainer.results.ka <- unlist(lapply(dagchainer.results.ka, function(x) ifelse(x > ka.max, NA, x)))
    dagchainer.results.ka <- unlist(lapply(dagchainer.results.ka, function(x) ifelse(x < ka.min, NA, x)))
    dagchainer.results.ks <- kaks$ks[match(dagchainer.pair.ids, kaks.pair.ids)]
    dagchainer.results.ks <- unlist(lapply(dagchainer.results.ks, function(x) ifelse(x > ks.max, NA, x)))
    dagchainer.results.ks <- unlist(lapply(dagchainer.results.ks, function(x) ifelse(x < ks.min, NA, x)))
    dagchainer.results.kaks <- dagchainer.results.ka / dagchainer.results.ks
    dagchainer.results.kaks[is.infinite(dagchainer.results.kaks)] <- NA
    dag <- dag %>% dplyr::mutate(ka = dagchainer.results.ka,
                                         ks = dagchainer.results.ks,
                                         kaks = dagchainer.results.kaks)
    dagchainer.results.split <- dag %>% group_by(gene1.chr, gene2.chr) %>% group_split(.keep = TRUE)
    g <- dag %>% group_by(gene1.chr, gene2.chr) %>% ggplot2::ggplot()
    if(colorBy == "Ka"){
      if(selfblast){
        g.ka <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = ka)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::geom_point(data = subset(g$data, gene1.chr == gene2.chr),
                              shape = 20, aes(x = gene1.mid, y = gene2.mid, col = ka)) +
          ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                      aes(slope = 1, intercept = 0)) +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.ka)
      }
      if(!selfblast){
        g.ka <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = ka)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.ka)
      }
    }
    if(colorBy == "Ks"){
      if(selfblast){
        g.ks <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = ks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::geom_point(data = subset(g$data, gene1.chr == gene2.chr),
                              shape = 20, aes(x = gene1.mid, y = gene2.mid, col = ks)) +
          ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                      aes(slope = 1, intercept = 0)) +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.ks)
      }
      if(!selfblast){
        g.ks <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid,
                                            col = ks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.ks)
      }
    }
    if(colorBy == "Ka/Ks"){
      if(selfblast){
        g.kaks <- g + ggplot2::geom_point(shape = 20,
                                          aes(x = gene2.mid, y = gene1.mid,
                                              col = kaks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::geom_point(data = subset(g$data, gene1.chr == gene2.chr),
                              shape = 20, aes(x = gene1.mid, y = gene2.mid, col = kaks)) +
          ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                      aes(slope = 1, intercept = 0)) +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.kaks)
      }
      if(!selfblast){
        g.kaks <- g + ggplot2::geom_point(shape = 20,
                                          aes(x = gene2.mid, y = gene1.mid,
                                              col = kaks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(DotPlotTitle)
        return(g.kaks)
      }
    }
  }
}
