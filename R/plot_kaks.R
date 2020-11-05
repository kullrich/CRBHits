#' @title plot_kaks
#' @name plot_kaks
#' @description This function plots Ka/Ks results obtained via `rbh2kaks()` function.
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [mandatory]
#' @param dag specify DAGchainer results as obtained via `rbh2dagchainer()` [default: NULL]
#' @param gene.position.cds1 specify gene position for cds1 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param gene.position.cds2 specify gene position for cds2 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param tandem.dups.cds1 specify tandem duplicates for cds1 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param tandem.dups.cds2 specify tandem duplicates for cds2 sequences (see \code{\link[CRBHits]{tandemdups}}) [default: NULL]
#' @param PlotTitle specify Plot title [default: Ka/Ks results]
#' @param PlotType specify Plot type: "h" histogram or "d" dotplot [default: h]
#' @param binw specify binwidth (see \code{\link[ggplot2]{geom_histogram}}) [default: 0.05]
#' @param splitByChr specify if plot should be split by chromosome [default: FALSE]
#' @param colorBy specify if Ka/Ks gene pairs should be colored by "rbh_class", dagchainer", "tandemdups" or "none" [default: rbh_class]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @param select.chr filter results for chromosome names [default: NULL]
#' @param doPlot specify plot [default: TRUE]
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap scale_colour_manual scale_colour_continuous aes geom_histogram ggtitle
#' @importFrom gridExtra grid.arrange
#' @examples
#' ## load example sequence data
#' data("ath_aly_ncbi_kaks", package="CRBHits")
#' ## plot Ka/Ks values - default
#' g <- plot_kaks(ath_aly_ncbi_kaks)
#' @export plot_kaks
#' @author Kristian K Ullrich

plot_kaks <- function(kaks,
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
                      select.chr = NULL,
                      doPlot = TRUE){
  gene1.chr <- NULL
  gene2.chr <- NULL
  gene1.mid <- NULL
  gene2.mid <- NULL
  ka <- NULL
  ks <- NULL
  rbh_class <- NULL
  if(attributes(kaks)$CRBHits.class != "kaks"){
    stop("Please obtain Ka/Ks via the rbh2kaks() function or add a 'kaks' class attribute")
  }
  selfblast <- attributes(kaks)$selfblast
  if(!is.null(dag)){
    if(attributes(dag)$CRBHits.class != "dagchainer"){
      stop("Please obtain DAGchainer results via the rbh2dagchainer() function or add a 'dagchainer' class attribute")
    }
    dagchainer.pair.ids <- apply(cbind(dag$gene1.seq.id, dag$gene2.seq.id),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    kaks.pair.ids <- apply(cbind(kaks$aa1, kaks$aa2),1, function(x) paste0(sort(x)[1],":",sort(x)[2]))
    kaks.results.dag <- dag$dagchainer_group[match(kaks.pair.ids, dagchainer.pair.ids)]
    kaks <- kaks %>% dplyr::mutate(dag = kaks.results.dag)
  }
  if(!is.null(gene.position.cds1)){
    if(attributes(gene.position.cds1)$CRBHits.class != "genepos"){
      stop("Please obtain gene position via the cds2genepos() function or add a 'genepos' class attribute")
    }
  }
  if(!is.null(gene.position.cds2)){
    if(attributes(gene.position.cds2)$CRBHits.class != "genepos"){
      stop("Please obtain gene position via the cds2genepos() function or add a 'genepos' class attribute")
    }
  }
  if(!is.null(gene.position.cds1) & !is.null(gene.position.cds2)){
    #add gene position
    gene.position.cds1.match <- gene.position.cds1[match(kaks$aa1,gene.position.cds1$gene.seq.id),]
    gene.position.cds2.match <- gene.position.cds2[match(kaks$aa2,gene.position.cds2$gene.seq.id),]
    kaks <- kaks %>% dplyr::mutate(gene1.seq.id = gene.position.cds1.match$gene.seq.id,
                                   gene1.chr = gene.position.cds1.match$gene.chr,
                                   gene1.start = gene.position.cds1.match$gene.start,
                                   gene1.end = gene.position.cds1.match$gene.end,
                                   gene1.mid = gene.position.cds1.match$gene.mid,
                                   gene1.strand = gene.position.cds1.match$gene.strand,
                                   gene1.idx = gene.position.cds1.match$gene.idx,
                                   gene2.seq.id = gene.position.cds2.match$gene.seq.id,
                                   gene2.chr = gene.position.cds2.match$gene.chr,
                                   gene2.start = gene.position.cds2.match$gene.start,
                                   gene2.end = gene.position.cds2.match$gene.end,
                                   gene2.mid = gene.position.cds2.match$gene.mid,
                                   gene2.strand = gene.position.cds2.match$gene.strand,
                                   gene2.idx = gene.position.cds2.match$gene.idx)
  }
  if(!is.null(tandem.dups.cds1)){
    if(attributes(tandem.dups.cds1)$CRBHits.class != "tandemdups"){
      stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
    }
  }
  if(!is.null(tandem.dups.cds2)){
    if(attributes(tandem.dups.cds2)$CRBHits.class != "tandemdups"){
      stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
    }
  }
  if(!is.null(tandem.dups.cds1) & !is.null(tandem.dups.cds2)){
    #add tandem dups
    tandem.dups.cds1.match <- tandem.dups.cds1[match(kaks$aa1,tandem.dups.cds1$gene.seq.id),]
    tandem.dups.cds2.match <- tandem.dups.cds2[match(kaks$aa2,tandem.dups.cds2$gene.seq.id),]
    kaks <- kaks %>% dplyr::mutate(gene1.tandem_group = tandem.dups.cds1.match$tandem_group,
                                   gene2.tandem_group = tandem.dups.cds2.match$tandem_group)
    tandem_group <- apply(cbind(kaks$gene1.tandem_group,kaks$gene2.tandem_group),1,function(x) paste0(x[1],":",x[2]))
    tandem_group[which(tandem_group == "NA:NA")] <- NA
    kaks <- kaks %>% dplyr::mutate(tandem_group = tandem_group)
  }
  if(!is.null(select.chr)){
    if(is.null(gene.position.cds1) | is.null(gene.position.cds2)){
      stop("Please provide gene position information. Can be obtained via the cds2genepos() function")
    }
    #filter for select.chr
    kaks <- kaks %>% dplyr::filter(gene1.chr %in% select.chr) %>% dplyr::filter(gene2.chr %in% select.chr)
  }
  if(PlotType == "d"){
    if(is.null(gene.position.cds1) | is.null(gene.position.cds2)){
      stop("Please provide gene position information. Can be obtained via the cds2genepos() function")
    }
  }
  if(splitByChr){
    if(is.null(gene.position.cds1) | is.null(gene.position.cds2)){
      stop("Please provide gene position information. Can be obtained via the cds2genepos() function")
    }
    kaks.split <- kaks %>% group_by(gene1.chr, gene2.chr) %>% group_split(.keep = TRUE)
  }
  if(colorBy == "dagchainer" | colorBy == "tandemdups"){
    if(is.null(gene.position.cds1) | is.null(gene.position.cds2)){
      stop("Please provide gene position information. Can be obtained via the cds2genepos() function")
    }
    if(colorBy == "dagchainer"){
      if(is.null(dag)){
        stop("Please obtain DAGchainer results via the rbh2dagchainer() function or add a 'dagchainer' class attribute")
      }
    }
    if(colorBy == "tandemdups"){
      if(is.null(tandem.dups.cds1)){
        stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
      }
      if(is.null(tandem.dups.cds2)){
        stop("Please obtain tandem duplicates via the tandemdups() function or add a 'tandemdups' class attribute")
      }
    }
  }
  #add kaks and filter for ka and ks min and max values
  if(attributes(kaks)$CRBHits.model == "Li"){
    kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x > ka.max, NA, x)))
    kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x < ka.min, NA, x)))
    kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x > ks.max, NA, x)))
    kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x < ks.min, NA, x)))
    kaks <- kaks %>% dplyr::mutate(kaks = kaks$ka / kaks$ks)
    kaks$kaks[is.infinite(kaks$kaks)] <- NA
  }
  if(attributes(kaks)$CRBHits.model == "YN"){
    suppressWarnings(kaks$ka <- as.numeric(kaks$ka))
    kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x > ka.max, NA, x)))
    kaks$ka <- unlist(lapply(kaks$ka, function(x) ifelse(x < ka.min, NA, x)))
    suppressWarnings(kaks$ks <- as.numeric(kaks$ks))
    kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x > ks.max, NA, x)))
    kaks$ks <- unlist(lapply(kaks$ks, function(x) ifelse(x < ks.min, NA, x)))
    kaks <- kaks %>% dplyr::mutate(kaks = kaks$ka / kaks$ks)
    kaks$kaks[is.infinite(kaks$kaks)] <- NA
  }
  if(colorBy == "none"){
    if(PlotType == "h"){
      if(!splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = kaks)) +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(PlotTitle)
        g.ka <- kaks %>% subset(!is.na(ka)) %>% ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka)) +
          ggplot2::ggtitle("Ka")
        g.ks <- kaks %>% subset(!is.na(ks)) %>% ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks)) +
          ggplot2::ggtitle("Ks")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
      if(splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = kaks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_continuous(type = "viridis") +
          ggplot2::ggtitle(PlotTitle)
        g.ka <- kaks %>% subset(!is.na(ka)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::ggtitle("Ka")
        g.ks <- kaks %>% subset(!is.na(ks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::ggtitle("Ks")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
    }
    if(PlotType == "d"){
      g <- kaks %>% subset(!is.na(ka)) %>%
        subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
        dplyr::group_by(gene1.chr, gene2.chr) %>%
        ggplot2::ggplot()
      g.kaks <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid, col = kaks)) +
        ggplot2::geom_point(shape = 20,
                            aes(x = gene1.mid, y = gene2.mid, col = kaks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_continuous(type = "viridis") +
        ggplot2::ggtitle(PlotTitle)
      g.ka <- g + ggplot2::geom_point(shape = 20,
                                        aes(x = gene2.mid, y = gene1.mid, col = ka)) +
        ggplot2::geom_point(shape = 20,
                            aes(x = gene1.mid, y = gene2.mid, col = ka)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_continuous(type = "viridis") +
        ggplot2::ggtitle("Ka")
      g.ks <- g + ggplot2::geom_point(shape = 20,
                                      aes(x = gene2.mid, y = gene1.mid, col = ks)) +
        ggplot2::geom_point(shape = 20,
                            aes(x = gene1.mid, y = gene2.mid, col = ks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_continuous(type = "viridis") +
        ggplot2::ggtitle("Ks")
      out <- list(g.kaks, g.ka, g.ks)
      names(out) <- c("g.kaks", "g.ka", "g.ks")
      if(doPlot){
        gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
      }
      return(out)
    }
  }
  if(colorBy == "rbh_class"){
    if(PlotType == "h"){
      if(!splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = rbh_class)) +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$rbh_class))), 50), na.value = "grey60") +
          ggplot2::ggtitle(paste0(PlotTitle, " - rbh_class"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>% ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = rbh_class)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$rbh_class))), na.value = "grey60") +
          ggplot2::ggtitle("Ka - rbh_class")
        g.ks <- kaks %>% subset(!is.na(ks)) %>% ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = rbh_class)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$rbh_class))), na.value = "grey60") +
          ggplot2::ggtitle("Ks - rbh_class")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
      if(splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = rbh_class)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$rbh_class))), 50), na.value = "grey60") +
          ggplot2::ggtitle(paste0(PlotTitle, "- rbh_class"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = rbh_class)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$rbh_class))), na.value = "grey60") +
          ggplot2::ggtitle("Ka - rbh_class")
        g.ks <- kaks %>% subset(!is.na(ks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = rbh_class)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$rbh_class))), na.value = "grey60") +
          ggplot2::ggtitle("Ks - rbh_class")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
    }
    if(PlotType == "d"){
      g <- kaks %>% subset(!is.na(ka)) %>%
        subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
        dplyr::group_by(gene1.chr, gene2.chr) %>%
        ggplot2::ggplot()
      g.kaks <- g + ggplot2::geom_point(shape = 21,
                                        aes(x = gene2.mid, y = gene1.mid, col = rbh_class, fill = kaks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = rbh_class, fill = kaks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$rbh_class))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::ggtitle(paste0(PlotTitle, " - rbh_class"))
      g.ka <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = rbh_class, fill = ka)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = rbh_class, fill = ka)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$rbh_class))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::ggtitle("Ka - rbh_class")
      g.ks <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = rbh_class, fill =ks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = rbh_class, fill = ks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$rbh_class))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::ggtitle("Ks - rbh_class")
      out <- list(g.kaks, g.ka, g.ks)
      names(out) <- c("g.kaks", "g.ka", "g.ks")
      if(doPlot){
        gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
      }
      return(out)
    }
  }
  if(colorBy == "dagchainer"){
    if(PlotType == "h"){
      if(!splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = dag)) +
          ggplot2::geom_point(data = subset(g$data, !is.na(dag)),
                              shape = 20, aes(x = ks, y = ka, col = dag)) +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$dag))), 50), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle(paste0(PlotTitle, " - dagchainer"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>% ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = dag)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ka - dagchainer")
        g.ks <- kaks %>% subset(!is.na(ks)) %>% ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = dag)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ks - dagchainer")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
      if(splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = dag)) +
          ggplot2::geom_point(data = subset(g$data, !is.na(dag)),
                              shape = 20, aes(x = ks, y = ka, col = dag)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$dag))), 50), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle(paste0(PlotTitle, " - dagchainer"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = dag)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ka - dagchainer")
        g.ks <- kaks %>% subset(!is.na(ks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = dag)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$dag))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ks - dagchainer")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
    }
    if(PlotType == "d"){
      g <- kaks %>% subset(!is.na(ka)) %>%
        subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
        dplyr::group_by(gene1.chr, gene2.chr) %>%
        ggplot2::ggplot()
      g.kaks <- g + ggplot2::geom_point(shape = 21,
                                        aes(x = gene2.mid, y = gene1.mid, col = dag, fill = kaks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = dag, fill = kaks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$dag))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle(paste0(PlotTitle, " - dagchainer"))
      g.ka <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = dag, fill = ka)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = dag, fill = ka)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$dag))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle("Ka - dagchainer")
      g.ks <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = dag, fill =ks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = dag, fill = ks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$dag))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle("Ks - dagchainer")
      out <- list(g.kaks, g.ka, g.ks)
      names(out) <- c("g.kaks", "g.ka", "g.ks")
      if(doPlot){
        gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
      }
      return(out)
    }
  }
  if(colorBy == "tandemdups"){
    if(PlotType == "h"){
      if(!splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = tandem_group)) +
          ggplot2::geom_point(data = subset(g$data, !is.na(tandem_group)),
                              shape = 20, aes(x = ks, y = ka, col = tandem_group)) +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$tandem_group))), 50), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle(paste0(PlotTitle, " - tandemdups"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>% ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = tandem_group)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$tandem_group))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ka - tandemdups")
        g.ks <- kaks %>% subset(!is.na(ks)) %>% ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = tandem_group)) +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$tandem_group))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ks - tandemdups")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
      if(splitByChr){
        g <- kaks %>% subset(!is.na(ka)) %>%
          subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.kaks <- g + ggplot2::geom_point(shape = 20, aes(x = ks, y = ka, col = tandem_group)) +
          ggplot2::geom_point(data = subset(g$data, !is.na(tandem_group)),
                              shape = 20, aes(x = ks, y = ka, col = tandem_group)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_colour_manual(
            values = col2transparent(CRBHitsColors(length(unique(kaks$tandem_group))), 50), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle(paste0(PlotTitle, " - tandemdups"))
        g.ka <- kaks %>% subset(!is.na(ka)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ka <- g.ka + ggplot2::geom_histogram(binwidth = binw, aes(x = ka, fill = tandem_group)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$tandem_group))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ka - tandemdups")
        g.ks <- kaks %>% subset(!is.na(ks)) %>%
          dplyr::group_by(gene1.chr, gene2.chr) %>%
          ggplot2::ggplot()
        g.ks <- g.ks + ggplot2::geom_histogram(binwidth = binw, aes(x = ks, fill = tandem_group)) +
          ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
          ggplot2::scale_fill_manual(
            values = CRBHitsColors(length(unique(kaks$tandem_group))), na.value = "grey60") +
          ggplot2::theme(legend.position="none") +
          ggplot2::ggtitle("Ks - tandemdups")
        out <- list(g.kaks, g.ka, g.ks)
        names(out) <- c("g.kaks", "g.ka", "g.ks")
        if(doPlot){
          gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
        }
        return(out)
      }
    }
    if(PlotType == "d"){
      g <- kaks %>% subset(!is.na(ka)) %>%
        subset(!is.na(ka)) %>% subset(!is.na(kaks)) %>%
        dplyr::group_by(gene1.chr, gene2.chr) %>%
        ggplot2::ggplot()
      g.kaks <- g + ggplot2::geom_point(shape = 21,
                                        aes(x = gene2.mid, y = gene1.mid, col = tandem_group, fill = kaks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = tandem_group, fill = kaks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$tandem_group))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle(paste0(PlotTitle, " - tandemdups"))
      g.ka <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = tandem_group, fill = ka)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = tandem_group, fill = ka)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$tandem_group))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle("Ka - tandemdups")
      g.ks <- g + ggplot2::geom_point(shape = 21,
                                      aes(x = gene2.mid, y = gene1.mid, col = tandem_group, fill =ks)) +
        ggplot2::geom_point(shape = 21,
                            aes(x = gene1.mid, y = gene2.mid, col = tandem_group, fill = ks)) +
        ggplot2::facet_grid(c("gene2.chr", "gene1.chr"), scales = "free") +
        ggplot2::geom_abline(data = subset(g$data, gene1.chr == gene2.chr),
                             aes(slope = 1, intercept = 0)) +
        ggplot2::scale_colour_manual(
          values = col2transparent(CRBHitsColors(length(unique(kaks$tandem_group))), 50), na.value = "grey60") +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle("Ks - tandemdups")
      out <- list(g.kaks, g.ka, g.ks)
      names(out) <- c("g.kaks", "g.ka", "g.ks")
      if(doPlot){
        gE <- gridExtra::grid.arrange(g.kaks, g.ka, g.ks, nrow = 1)
      }
      return(out)
    }
  }
}
