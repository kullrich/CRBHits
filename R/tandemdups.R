#' @title tandemdups
#' @name tandemdups
#' @description This function assigns tandem duplicates given (conditional-)reciprocal best hit (CRBHit) pairs and their chromosomal gene position.
#' The function is ported into R from Haug-Baltzell et al. (2017) https://github.com/LyonsLab/coge/blob/master/bin/dagchainer/tandems.py
#' @param rbhpairs (conditional-)reciprocal best hit (CRBHit) pair result (see \code{\link[CRBHits]{cds2rbh}}) [mandatory]
#' @param genepos Gene position matrix as obtained from \code{cds2genepos} [mandatory]
#' @param dupdist Maximum distance between 2 gene positions on the same chromosome which will call them as a pair of local duplicates. If they are farther apart than \code{dipdist}, but share a common hit within \code{dupdist} of both, then they will still be in the same set of local duplicates. [default: 5]
#' @return \code{matrix}
#' 1: $gene.seq.id\cr
#' 2: $gene.chr\cr
#' 3: $gene.start\cr
#' 4: $gene.end\cr
#' 5: $gene.mid\cr
#' 6: $gene.strand\cr
#' 7: $gene.idx\cr
#' 8: $tandem_group\cr
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys
#' @importFrom foreach foreach %do%
#' @importFrom stringr str_sub
#' @seealso \code{\link[CRBHits]{isoform2longest}}
#' @references Haug-Beltzell A et al. (2017) SynMap2 and SynMap3D: web-based whole-genome synteny browsers. \emph{Bioinformatics} \bold{33(14)}, 2197-2198.
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' ## get selfhits CRBHit pairs
#' ath_selfhits_crbh <- cds2rbh(ath, ath, plotCurve = TRUE)
#' ## get gene position
#' ath.genepos <- cds2genepos(ath, "ENSEMBL")
#' ## get tandem duplicate results
#' ath_selfblast_crbh.tandemdups <- tandemdups(ath_selfhits_crbh,
#'                                             ath.genepos)
#' head(ath_selfblast_crbh.tandemdups)
#' @export tandemdups
#' @author Kristian K Ullrich
tandemdups <- function(rbhpairs, genepos, dupdist = 5){
  if(attributes(rbhpairs)$CRBHits.class != "crbh"){
    stop("Please obtain rbhpairs via the cds2rbh() or the cdsfile2rbh() function")
  }
  if(attributes(genepos)$CRBHits.class != "genepos"){
    stop("Please obtain gene position via the cds2genepos() function or add a genepos class attribute")
  }
  if(dupdist < 1){
    stop("dupdist must be at least 1")
  }
  localdups <- function(selfhits_per_chr, dupdist){
    #define tandem duplicate group prefix
    tandem.group.prefix <- paste0("TDG-",unique(selfhits_per_chr$gene1.chr))
    #parent duplicate is defined as the first hit based on gene.idx
    #if gene1.idx == gene2.idx continue
    selfhits_per_chr <- selfhits_per_chr %>% dplyr::arrange(gene1.idx, gene2.idx)
    selfhits_per_chr.distinct <- selfhits_per_chr %>% dplyr::distinct(gene1.idx, gene2.idx, .keep_all = TRUE)
    #get gene distance
    selfhits_per_chr.distinct <- selfhits_per_chr.distinct %>% dplyr::mutate(
      gene.dist = abs(selfhits_per_chr$gene1.idx - selfhits_per_chr$gene2.idx))
    #retain only hits with smaller gene distance as dupdist
    selfhits_per_chr.distinct.dupdist <- selfhits_per_chr.distinct %>% dplyr::filter(gene.dist <= dupdist)
    #if sort(gene1.idx, gene2.idx) is duplicated continue
    selfhits_per_chr.distinct.dupdist <- selfhits_per_chr.distinct.dupdist %>%
      dplyr::mutate(seendups = apply(
        cbind(selfhits_per_chr.distinct.dupdist$gene1.idx, selfhits_per_chr.distinct.dupdist$gene2.idx), 1,
        function(x) paste0(sort(x),collapse = ":"))) %>% dplyr::distinct(seendups, .keep_all = TRUE)
    #define local duplicate groups
    loc_dups <- list()
    inserted <- list()
    out <- foreach::foreach(qloc = selfhits_per_chr.distinct.dupdist$gene1.idx,
                     sloc = selfhits_per_chr.distinct.dupdist$gene2.idx, .combine=c) %do% {
      qloc <- as.character(qloc)
      sloc <- as.character(sloc)
      qloc_sloc <- paste0(sort(c(qloc, sloc)), collapse = ":")
      qloc_sloc_used <- FALSE
      if(qloc %in% names(inserted)){
        loc_dups[[inserted[[qloc]]]] <- c(loc_dups[[inserted[[qloc]]]], sloc)
        inserted[[sloc]] <- inserted[[qloc]]
        qloc_sloc_used <- TRUE
      } else if(qloc %in% names(loc_dups)){
        loc_dups[[qloc]] <- c(loc_dups[[qloc]], sloc)
        inserted[[sloc]] <- qloc
        qloc_sloc_used <- TRUE
      }
      if(sloc %in% names(inserted)){
        loc_dups[[inserted[[sloc]]]] <- c(loc_dups[[inserted[[sloc]]]], qloc)
        inserted[[qloc]] <- inserted[[sloc]]
        qloc_sloc_used <- TRUE
      } else if(sloc %in% names(loc_dups)){
        loc_dups[[sloc]] <- c(loc_dups[[sloc]], qloc)
        inserted[[qloc]] <- sloc
        qloc_sloc_used <- TRUE
      }
      if(!qloc_sloc_used){
        loc_dups[[qloc]] <- sloc
      }
      return(NULL)}
    #get loc_dup names as parent duplicate
    loc_dups.names <- names(loc_dups)[order(as.numeric(names(loc_dups)))]
    #add parent duplicate into loc_dups groups
    for(lkey in loc_dups.names){loc_dups[[lkey]] <- sort(unique(as.numeric(c(loc_dups[[lkey]], lkey))))}
    kidx1 <- 0
    kidx2 <- 1
    while( (kidx1 + kidx2) < length(loc_dups.names) ){
      k1 <- loc_dups.names[kidx1 + 1]
      vals1 <- as.numeric(unique(loc_dups[[k1]]))
      #so any other dupset within dupdist of this one could
      #potentially have a common dup, so look through all of them.
      while ( kidx2 < dupdist & (kidx1 + kidx2) < length(loc_dups.names) ){
        k2 <- loc_dups.names[kidx1 + kidx2 + 1]
        vals2 <- as.numeric(unique(loc_dups[[k2]]))
        # they dont have any in common.
        if(!any(vals1 %in% vals2)){
          kidx2 <- kidx2 + 1
          next()
        }
        #can do destructive inside loop because we delete k1 not k2
        #and we're using a separate copy of the keys.
        if(any(vals1 %in% vals2)){
          loc_dups[[k1]] <- NULL
          loc_dups[[k2]] <- as.numeric(unique(c(loc_dups[[k2]], vals1)))
          break()
        }
      }
      kidx1 <- kidx1 + 1
      kidx2 <- 1
    }
    #return loc_dups
    #remove repeats and parents and sort.
    #use new dups to make sure lowest number accn is the parent.
    if(length(loc_dups) == 0){
      out <- selfhits_per_chr %>% dplyr::mutate(
        gene1.tandem_group = NA,
        gene2.tandem_group = NA)
      return(out)
    } else{
      tandem.groups <- loc_dups
      tandem.groups <- lapply(tandem.groups, unique)
      names(tandem.groups) <- as.character(unlist(lapply(tandem.groups, min)))
      tandem.groups <- tandem.groups[order(as.numeric(names(tandem.groups)))]
      names(tandem.groups) <- paste0(tandem.group.prefix, "-", as.character(unlist(lapply(tandem.groups, min))))
      tandem.groups.stack <- stack(tandem.groups)
      tandem.groups <- setNames(tandem.groups.stack$values, tandem.groups.stack$ind)
      out <- selfhits_per_chr %>% dplyr::mutate(
        gene1.tandem_group = names(tandem.groups[match(selfhits_per_chr$gene1.idx, tandem.groups)]),
        gene2.tandem_group = names(tandem.groups[match(selfhits_per_chr$gene2.idx, tandem.groups)]))
    return(out)
    }
  }
  matchList <- dplyr::bind_cols(
    genepos[match(rbhpairs$crbh.pairs$aa1,
                  genepos$gene.seq.id), , drop = FALSE] %>% select(gene1.seq.id = gene.seq.id,
                                                                   gene1.chr = gene.chr,
                                                                   gene1.start = gene.start,
                                                                   gene1.end = gene.end,
                                                                   gene1.mid = gene.mid,
                                                                   gene1.strand = gene.strand,
                                                                   gene1.idx = gene.idx),
    genepos[match(rbhpairs$crbh.pairs$aa2,
                  genepos$gene.seq.id), , drop = FALSE] %>% select(gene2.seq.id = gene.seq.id,
                                                                   gene2.chr = gene.chr,
                                                                   gene2.start = gene.start,
                                                                   gene2.end = gene.end,
                                                                   gene2.mid = gene.mid,
                                                                   gene2.strand = gene.strand,
                                                                   gene2.idx = gene.idx))
  #group by chromosome and CRBHit pairs on the same chromosome
  selfhits <- matchList %>% dplyr::group_by(gene1.chr) %>%
    dplyr::filter(gene1.chr == gene2.chr) %>% dplyr::arrange(gene1.chr, gene1.idx, gene2.idx)
  #get loc_dups by group
  selfhits.loc_dups <- selfhits %>% dplyr::group_map(~ localdups(.x, dupdist), .keep = TRUE)
  names(selfhits.loc_dups) <- unlist(dplyr::group_keys(selfhits))
  tandemdups <- do.call(rbind, selfhits.loc_dups)
  tandemdups <- rbind(tandemdups %>% dplyr::select(gene.seq.id = gene1.seq.id,
   gene.chr = gene1.chr, gene.start = gene1.start, gene.end = gene1.end,
   gene.mid = gene1.mid, gene.strand = gene1.strand, gene.idx = gene1.idx,
   tandem_group = gene1.tandem_group),
   tandemdups %>% dplyr::select(gene.seq.id = gene2.seq.id,
   gene.chr = gene2.chr, gene.start = gene2.start, gene.end = gene2.end,
   gene.mid = gene2.mid, gene.strand = gene2.strand, gene.idx = gene2.idx,
   tandem_group = gene2.tandem_group)) %>% dplyr::distinct(gene.seq.id, .keep_all =TRUE) %>%
   dplyr::filter(!is.na(tandem_group)) %>% dplyr::arrange(gene.chr, gene.idx)
  attr(tandemdups, "CRBHits.class") <- "tandemdups"
  return(tandemdups)
}
