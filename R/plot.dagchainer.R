#' @title plot.dagchainer
#' @name plot.dagchainer
#' @description This function plots DAGchainer (http://dagchainer.sourceforge.net/) results obtained via `rbh2dagchainer()` function.
#' @param dag specify DAGchainer results as obtained via `rbh2dagchainer()` [mandatory]
#' @param DotPlotTitle specify DotPlot title [default: DAGchainer results]
#' @param colorBy specify if dagchainer groups should be colored by "Ka", "Ks", "Ka/KS" or "none" [default: none]
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [default: NULL]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap
#' @export plot.dagchainer
#' @author Kristian K Ullrich

plot.dagchainer <- function(dag, DotPlotTitle = "DAGchainer results", colorBy = "none", kaks = NULL, ka.max = 5, ks.max = 5, ka.min = 0, ks.min = 0){
  if(colorBy == "none"){
    dagchainer.results.split <- dagchainer.results %>% group_by(gene1.chr, gene2.chr) %>% group_split(.keep = TRUE)
    g <- dagchainer.results %>% group_by(gene1.chr, gene2.chr) %>% ggplot2::ggplot()
    g1 <- g + ggplot2::geom_point(shape = 20,
                                  aes(x = gene2.mid, y = gene1.mid,
                                      col = as.factor(dagchainer_group))) +
      ggplot2::theme(legend.position="none") +
      ggplot2::facet_grid(c("gene1.chr", "gene2.chr"), scales = "free") +
      ggplot2::scale_colour_manual(
        values = CRBHitsColors(length(unique(dagchainer.results$dagchainer_group))))
    g2 <- g1 + ggplot2::geom_point(data = subset(g1$data, gene1.chr == gene2.chr), shape = 20,
                                   aes(x = gene1.mid, y = gene2.mid,
                                       col = as.factor(dagchainer_group))) +
      geom_abline(data = subset(g1$data, gene1.chr == gene2.chr), aes(slope = 1, intercept = 0))
    g2 + ggplot2::ggtitle(DotPlotTitle)
  }
  if(colorBy == "Ka" | colorBy == "Ks" | colorBy == "Ka/Ks"){
    dagchainer.pair.ids <- apply(cbind(dagchainer.results$gene1.seq.id, dagchainer.results$gene2.seq.id),1, function(x) paste0(x[1],":",x[2]))
    kaks.pair.ids <- apply(cbind(kaks$aa1, kaks$aa2),1, function(x) paste0(x[1],":",x[2]))

        dagchainer.results.gene1.ka <- kaks$ka[match(dagchainer.results$gene1.seq.id, kaks$aa1)]
    dagchainer.results.gene2.ka <- kaks$ka[match(dagchainer.results$gene2.seq.id, kaks$aa2)]
  }
}
