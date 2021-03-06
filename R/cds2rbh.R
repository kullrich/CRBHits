#' @title cds2rbh
#' @name cds2rbh
#' @description This function calculates (conditional-)reciprocal best hit (CRBHit) pairs from two CDS \code{DNAStringSet}'s.
#' CRBHit pairs were introduced by \emph{Aubry S, Kelly S et al. (2014)}.
#' Sequence searches are performed with \bold{last} \emph{Kiełbasa, SM et al. (2011)}.
#' If one specifies cds1 and cds2 as the same input a selfblast is conducted.
#' @param cds1 cds1 sequences as \code{DNAStringSet} [mandatory]
#' @param cds2 cds2 sequences as \code{DNAStringSet} [mandatory]
#' @param lastpath specify the PATH to the last binaries [default: /extdata/last-1243/bin/]
#' @param lastD last option D: query letters per random alignment [default: 1e6]
#' @param outpath specify the output PATH [default: /tmp]
#' @param crbh specify if conditional-reciprocal hit pairs should be retained as secondary hits [default: TRUE]
#' @param keepSingleDirection specify if single direction secondary hit pairs should be retained [default: FALSE]
#' @param eval evalue [default: 1e-3]
#' @param qcov query coverage [default: 0.0]
#' @param tcov target coverage [default: 0.0]
#' @param pident percent identity [default: 0.0]
#' @param alnlen alignment length [default: 0.0]
#' @param rost1999 specify if hit pairs should be filter by equation 2 of Rost 1999 [default: FALSE]
#' @param filter specify additional custom filters as list to be applied on hit pairs [default: NULL]
#' @param plotCurve specify if crbh fitting curve should be plotted [default: FALSE]
#' @param fit.type specify if mean or median should be used for fitting [default: mean]
#' @param fit.varweight factor for fitting function to consider neighborhood [default: 0.1]
#' @param fit.min specify minimum neighborhood alignment length [default: 5]
#' @param longest.isoform specify if cds sequences should be removed to the longest isoform (only possible if data was accessed from NCBI or ENSEMBL) [default: FALSE]
#' @param isoform.source specify cds sequences source (either NCBI or ENSEMBL) [default: NCBI]
#' @param threads number of parallel threads [default: 1]
#' @param remove specify if last result files should be removed [default: TRUE]
#' @return List of three (crbh = FALSE)\cr
#' 1: $crbh.pairs\cr
#' 2: $crbh1 matrix; query > target\cr
#' 3: $crbh2 matrix; target > query\cr
#' \cr
#' List of four (crbh = TRUE)\cr
#' 1: $crbh.pairs\cr
#' 2: $crbh1 matrix; query > target\cr
#' 3: $crbh2 matrix; target > query\cr
#' 4: $rbh1_rbh2_fit; evalue fitting function
#' @importFrom Biostrings writeXStringSet
#' @importFrom graphics legend par points
#' @importFrom stats median splinefun
#' @importFrom utils read.table
#' @importFrom tidyr %>%
#' @seealso \code{\link[CRBHits]{cdsfile2rbh}},
#' \code{\link[CRBHits]{isoform2longest}}
#' @references Aubry S, Kelly S et al. (2014) Deep Evolutionary Comparison of Gene Expression Identifies Parallel Recruitment of Trans-Factors in Two Independent Origins of C4 Photosynthesis. \emph{PLOS Genetics}, \bold{10(6)} e1004365.
#' @references Kiełbasa, SM et al. (2011) Adaptive seeds tame genomic sequence comparison. \emph{Genome research}, \bold{21(3)}, 487-493.
#' @references Rost B. (1999). Twilight zone of protein sequence alignments. \emph{Protein Engineering}, \bold{12(2)}, 85-94.
#' @examples
#' ## compile last-1243 within CRBHits
#' make_last()
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' data("aly", package="CRBHits")
#' ## get CRBHit pairs
#' ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE)
#' dim(ath_aly_crbh$crbh.pairs)
#' ## get classical reciprocal best hit (RBHit) pairs
#' ath_aly_rbh <- cds2rbh(ath, aly, crbh = FALSE)
#' dim(ath_aly_rbh$crbh.pairs)
#' ## filter for evalue 1e-100
#' ath_aly_crbh.eval100 <- cds2rbh(ath, aly, eval = 1e-100)
#' dim(ath_aly_crbh.eval100$crbh.pairs)
#' ## filter for query coverage
#' ath_aly_crbh.qcov <- cds2rbh(ath, aly, qcov = 0.5)
#' dim(ath_aly_crbh.qcov$crbh.pairs)
#' ## custom filter for e.g. bit score (column 12)
#' ## Note: multiple filters can be given in filter param as list
#' myfilter <- function(rbh, value = 500.0){
#' return(rbh[as.numeric(rbh[, 12]) >= value , , drop = FALSE])
#' }
#' ath_aly_crbh.custom <- cds2rbh(ath, aly, filter = list(myfilter))
#' dim(ath_aly_crbh.custom$crbh.pairs)
#' ## selfblast
#' ath_selfblast_crbh <- cds2rbh(ath, ath, plotCurve = TRUE)
#' @export cds2rbh
#' @author Kristian K Ullrich

cds2rbh <- function(cds1, cds2,
                    lastpath = paste0(find.package("CRBHits"),
                                              "/extdata/last-1243/bin/"),
                    lastD = 1e6,
                    outpath = "/tmp",
                    crbh = TRUE,
                    keepSingleDirection = FALSE,
                    eval = 1e-3,
                    qcov = 0.0,
                    tcov = 0.0,
                    pident = 0.0,
                    alnlen = 0.0,
                    rost1999 = FALSE,
                    filter = NULL,
                    plotCurve = FALSE,
                    fit.type = "mean",
                    fit.varweight = 0.1,
                    fit.min = 5,
                    longest.isoform = FALSE,
                    isoform.source = "NCBI",
                    threads = 1,
                    remove = TRUE){
  #internal function to fit evalue by length
  fitSpline <- function(alnlength, evalue, fit.type, fit.varweight, fit.min){
    log10evalue <- -log10(evalue)
    log10evalue[is.infinite(log10evalue)] <- 324
    x <- cbind(alnlength, log10evalue)
    x <- x[order(x[, 1]), ]
    x.max <- max(x[, 1])
    fitMatrix <- matrix(0, ncol = 2, nrow = x.max)
    for(i in seq(from = 1, to = x.max)){
      fitMatrix[i, 1] <- i
      s <- round(i * fit.varweight)
      if(s < fit.min){s <- fit.min}
      smin <- i - s
      smax <- i + s
      s.idx <- which(x[, 1] >= smin & x[, 1] <= smax)
      if(length(s.idx) == 0){s.value <- 0}
      if(length(s.idx) != 0){
        if(fit.type == "mean"){
          s.value <- mean(x[s.idx, 2])
        }
        if(fit.type == "median"){
          s.value <- median(x[s.idx, 2])
        }
      }
      if(i == 1){
        fitMatrix[i, 2] <- s.value
      }
      if(i > 1){
        if(fitMatrix[i - 1, 2] <= s.value){
          fitMatrix[i, 2] <- s.value
        }
        if(fitMatrix[i - 1, 2] > s.value){
          fitMatrix[i, 2] <- fitMatrix[i - 1, 2]
        }
      }
    }
    fitMatrixfun <- splinefun(fitMatrix[,1], fitMatrix[, 2])
    return(fitMatrixfun)
  }
  if(!dir.exists(lastpath)){stop("Error: last PATH does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_last() function.")}
  if(!file.exists(paste0(lastpath, "lastdb"))){stop("Error: lastdb binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_last() function.")}
  if(!file.exists(paste0(lastpath, "lastal"))){stop("Error: lastal binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make_last() function.")}
  selfblast <- FALSE
  #suppressWarnings since if length of cds1 and cds2 differ, this will raise a warning
  if(suppressWarnings(any(cds1 == cds2))){
    selfblast <- TRUE
  }
  aa1file <- tempfile("aa1_", outpath)
  aa2file <- tempfile("aa2_", outpath)
  aa1dbfile <- tempfile("aa1db_", outpath)
  aa2dbfile <- tempfile("aa2db_", outpath)
  aa2_aa1_lastout <- tempfile("aa2_aa1_lastout_", outpath)
  aa1_aa2_lastout <- tempfile("aa1_aa2_lastout_", outpath)
  if(longest.isoform){
    Biostrings::writeXStringSet(cds2aa(isoform2longest(cds1, isoform.source)), file = aa1file)
    Biostrings::writeXStringSet(cds2aa(isoform2longest(cds2, isoform.source)), file = aa2file)
  }
  if(!longest.isoform){
    Biostrings::writeXStringSet(cds2aa(cds1), file = aa1file)
    Biostrings::writeXStringSet(cds2aa(cds2), file = aa2file)
  }
  system(paste0(lastpath, "lastdb -p -cR01 -P ", threads," ", aa1dbfile, " ", aa1file))
  system(paste0(lastpath, "lastdb -p -cR01 -P ", threads," ", aa2dbfile, " ", aa2file))
  system(paste0(lastpath, "lastal -f BlastTab+ -P ", threads, " -D", lastD, " ", aa1dbfile, " ", aa2file, " > ", aa2_aa1_lastout))
  system(paste0(lastpath, "lastal -f BlastTab+ -P ", threads, " -D", lastD, " ", aa2dbfile, " ", aa1file, " > ", aa1_aa2_lastout))
  aa1_aa2 <- read.table(aa1_aa2_lastout, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  aa2_aa1 <- read.table(aa2_aa1_lastout, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(aa1_aa2) <- colnames(aa2_aa1) <- c("query_id", "subject_id", "perc_identity",
                                              "alignment_length", "mismatches",
                                              "gap_opens", "q_start", "q_end",
                                              "s_start", "s_end", "evalue",
                                              "bit_score", "query_length",
                                              "subject_length", "raw_score")
  if(remove){
    system(paste0("rm ", aa1file))
    system(paste0("rm ", aa2file))
    system(paste0("rm ", aa1dbfile, "*"))
    system(paste0("rm ", aa2dbfile, "*"))
    system(paste0("rm ", aa2_aa1_lastout))
    system(paste0("rm ", aa1_aa2_lastout))
  }
  #selfblast
  if(selfblast){
    aa1_aa2 <- aa1_aa2[aa1_aa2[, 1] != aa1_aa2[, 2], , drop = FALSE]
    aa2_aa1 <- aa2_aa1[aa2_aa1[, 1] != aa2_aa1[, 2], , drop = FALSE]
    if(dim(aa1_aa2)[1] == 0 & dim(aa2_aa1)[1] == 0){
      {stop("No recirpocal best hits!")}
    }
  }
  #apply standard filters on hit pairs
  aa1_aa2 <- aa1_aa2 %>% filter_eval(eval) %>% filter_qcov(qcov) %>%
                         filter_tcov(tcov) %>% filter_pident(pident) %>%
                         filter_alnlen(alnlen)
  aa2_aa1 <- aa2_aa1 %>% filter_eval(eval) %>% filter_qcov(qcov) %>%
                         filter_tcov(tcov) %>% filter_pident(pident) %>%
                         filter_alnlen(alnlen)
  if(rost1999){
    aa1_aa2 <- aa1_aa2 %>% filter_rost1999
    aa2_aa1 <- aa2_aa1 %>% filter_rost1999
  }
  #apply additional filters on hit pairs
  for(f_ in filter){
    aa1_aa2 <- aa1_aa2 %>% f_
    aa2_aa1 <- aa2_aa1 %>% f_
  }
  if(dim(aa1_aa2)[1] == 0 & dim(aa2_aa1)[1] == 0){
    {stop("No recirpocal best hits!")}
  }
  aa1_aa2.idx <- paste0(aa1_aa2[, 1], "\t", aa1_aa2[, 2])
  aa2_aa1.idx <- paste0(aa2_aa1[, 2], "\t", aa2_aa1[, 1])
  #deduplicate hit pairs and only retain the best hit per query
  aa1_aa2.dedup <- aa1_aa2[!duplicated(aa1_aa2[, 1]), , drop = FALSE]
  aa2_aa1.dedup <- aa2_aa1[!duplicated(aa2_aa1[, 1]), , drop = FALSE]
  aa1_aa2.dedup.idx <- paste0(aa1_aa2.dedup[, 1], "\t", aa1_aa2.dedup[, 2])
  aa2_aa1.dedup.idx <- paste0(aa2_aa1.dedup[, 2], "\t", aa2_aa1.dedup[, 1])
  #reduce to reciprocal best hits
  rbh1 <- aa1_aa2.dedup[which(aa1_aa2.dedup.idx %in% aa2_aa1.dedup.idx), , drop = FALSE]
  rbh2 <- aa2_aa1.dedup[which(aa2_aa1.dedup.idx %in% aa1_aa2.dedup.idx), , drop = FALSE]
  if(selfblast){
    rbh1 <- rbh1[!duplicated(apply(rbh1[, 1:2], 1, function(x) paste0(sort(x), collapse="\t"))), , drop = FALSE]
    rbh2 <- rbh2[!duplicated(apply(rbh2[, 1:2], 1, function(x) paste0(sort(x), collapse="\t"))), , drop = FALSE]
  }
  #if no crbh - done
  if(!crbh){
    rbh <- rbh1[, 1:2]
    rbh <- cbind(rbh, "rbh")
    colnames(rbh) <- c("aa1", "aa2", "rbh_class")
    out <- list(rbh, cbind(rbh1, "rbh"), cbind(rbh2, "rbh"))
    names(out) <- c("crbh.pairs", "crbh1", "crbh2")
    if(!selfblast){
      attr(out, "CRBHits.class") <- "crbh"
      attr(out, "crbh") <- FALSE
      attr(out, "keepSingleDirection") <- FALSE
      attr(out, "selfblast") <- selfblast
      return(out)
    }
    if(selfblast){
      attr(out, "CRBHits.class") <- "crbh"
      attr(out, "crbh") <- FALSE
      attr(out, "keepSingleDirection") <- FALSE
      attr(out, "selfblast") <- selfblast
      return(out)
    }
  }
  #if crbh - continue
  if(crbh){
    #fit evalue by length
    if(!selfblast){
      rbh1_rbh2_fit <- fitSpline(c(rbh1[, 4], rbh2[, 4]),
                                 c(rbh1[, 11], rbh2[, 11]),
                                 fit.type,
                                 fit.varweight,
                                 fit.min)
    }
    #selfblast
    if(selfblast){
      rbh1_rbh2_fit <- fitSpline(c(rbh1[, 4]),
                                 c(rbh1[, 11]),
                                 fit.type,
                                 fit.varweight,
                                 fit.min)
    }
    #internal function to filter hit pairs using rbh1_rbh2_fit
    filter_crbh <- function(x){
      minuslog10evalue_by_fit <- lapply(as.numeric(x[, 4]), rbh1_rbh2_fit)
      return(x[as.numeric(x[,16]) >= minuslog10evalue_by_fit, , drop = FALSE])
    }
    #remove reciprocal best hits from hit pairs to only look into secondary hits
    if(!selfblast){
      rbh1.idx <- paste0(rbh1[, 1], "\t", rbh1[, 2])
      rbh2.idx <- paste0(rbh2[, 2], "\t", rbh2[, 1])
    }
    if(selfblast){
      rbh1.idx <- c(paste0(rbh1[, 1], "\t", rbh1[, 2]), paste0(rbh1[, 2], "\t", rbh1[, 1]))
      rbh2.idx <- c(paste0(rbh2[, 2], "\t", rbh2[, 1]), paste0(rbh2[, 1], "\t", rbh2[, 2]))
    }
    aa1_aa2.red <- aa1_aa2[!aa1_aa2.idx %in% rbh1.idx, , drop = FALSE]
    aa2_aa1.red <- aa2_aa1[!aa2_aa1.idx %in% rbh2.idx, , drop = FALSE]
    #add -log10(evalue)
    aa1_aa2.red <- cbind(aa1_aa2.red, -log10(aa1_aa2.red[, 11]))
    aa1_aa2.red[is.infinite(aa1_aa2.red[, 16]), 16] <- 324
    aa2_aa1.red <- cbind(aa2_aa1.red, -log10(aa2_aa1.red[, 11]))
    aa2_aa1.red[is.infinite(aa2_aa1.red[, 16]), 16] <- 324
    #filter retained hit pairs with rbh1_rbh2_fit
    aa1_aa2.red <- filter_crbh(aa1_aa2.red)
    aa2_aa1.red <- filter_crbh(aa2_aa1.red)
    aa1_aa2.red.idx <- paste0(aa1_aa2.red[, 1], "\t", aa1_aa2.red[, 2])
    aa2_aa1.red.idx <- paste0(aa2_aa1.red[, 2], "\t", aa2_aa1.red[, 1])
    #deduplicate hit pairs and only retain the best hit per HSP
    aa1_aa2.red.dedup <- aa1_aa2.red[!duplicated(aa1_aa2.red.idx), , drop = FALSE]
    aa2_aa1.red.dedup <- aa2_aa1.red[!duplicated(aa2_aa1.red.idx), , drop = FALSE]
    #split into reciprocal direction secondary hits (rbh.sec) and single direction (single)
    aa1_aa2.red.dedup.idx <- paste0(aa1_aa2.red.dedup[, 1], "\t", aa1_aa2.red.dedup[, 2])
    aa2_aa1.red.dedup.idx <- paste0(aa2_aa1.red.dedup[, 2], "\t", aa2_aa1.red.dedup[, 1])
    rbh1.sec <- aa1_aa2.red.dedup[which(aa1_aa2.red.dedup.idx %in% aa2_aa1.red.dedup.idx), , drop = FALSE]
    rbh2.sec <- aa2_aa1.red.dedup[which(aa2_aa1.red.dedup.idx %in% aa1_aa2.red.dedup.idx), , drop = FALSE]
    if(selfblast){
      rbh1.sec <- rbh1.sec[!duplicated(apply(rbh1.sec[, 1:2], 1, function(x) paste0(sort(x), collapse="\t"))), , drop = FALSE]
      rbh2.sec <- rbh2.sec[!duplicated(apply(rbh2.sec[, 1:2], 1, function(x) paste0(sort(x), collapse="\t"))), , drop = FALSE]
    }
    single1 <- aa1_aa2.red.dedup[which(!aa1_aa2.red.dedup.idx %in% aa2_aa1.red.dedup.idx), , drop = FALSE]
    single2 <- aa2_aa1.red.dedup[which(!aa2_aa1.red.dedup.idx %in% aa1_aa2.red.dedup.idx), , drop = FALSE]
    #if plotCurve plot fitting
    if(plotCurve){
      len <- rbh1[, 4]
      log10alnlen <- log10(len)
      minuslog10evalue <- -log10(rbh1[, 11])
      minuslog10evalue[is.infinite(minuslog10evalue)] <- 324
      plot(x = log10alnlen, y = minuslog10evalue,
           pch = 20,
           main = "Accept / Reject secondary hits as homologs",
           ylab = "-log10(evalue)",
           xlab = "log10(alnlength)",
           col = col2transparent("#8EBCB5", 50),
           cex = 0.75)
      points(x = log10(rbh1.sec[, 4]), y = -log10(rbh1.sec[, 11]),
             pch = 21,
             bg = col2transparent("#4D83AB", 25),
             cex = 1)
      points(x = log10(single1[, 4]), y = -log10(single1[, 11]),
             pch = 21,
             bg = col2transparent("#CBC106", 25),
             cex = 1)
      if(selfblast){
        points(x = log10(1:max(rbh1[, 4])),
               y = rbh1_rbh2_fit(seq(from = 1, to = max(rbh1[, 4]))),
               type = "l",
               lwd = 2,
               col = "#9E163C")
        legend("bottomright",
               legend = c("rbh", "sec", "single"),
               col = c("#8EBCB5",
                       col2transparent("#4D83AB", 25),
                       col2transparent("#CBC106", 25)),
               bty = "n",
               pch = 20)
      }
      if(!selfblast){
        points(x = log10(single2[, 4]), y = -log10(single2[, 11]),
               pch = 21,
               bg = col2transparent("#CB7B26", 25),
               cex = 1)
        points(x = log10(1:max(rbh1[, 4])),
               y = rbh1_rbh2_fit(seq(from = 1, to = max(rbh1[, 4]))),
               type = "l",
               lwd = 2,
               col = "#9E163C")
        legend("bottomright",
               legend = c("rbh", "sec", "single1", "single2"),
               col = c("#8EBCB5",
                       col2transparent("#4D83AB", 25),
                       col2transparent("#CBC106", 25),
                       col2transparent("#CB7B26", 25)),
               bty = "n",
               pch = 20)
      }
    }
    #if no keepSingleDirection - done
    if(!keepSingleDirection){
      crbh1 <- data.frame(Map(c ,cbind(rbh1, "rbh"), cbind(rbh1.sec[, 1:15], "sec")))
      colnames(crbh1)[16] <- "rbh_class"
      crbh2 <- data.frame(Map(c ,cbind(rbh2, "rbh"), cbind(rbh2.sec[, 1:15], "sec")))
      colnames(crbh2)[16] <- "rbh_class"
      crbh <- crbh1[, c(1:2,16)]
      colnames(crbh) <- c("aa1", "aa2", "rbh_class")
      out <- list(crbh, crbh1, crbh2, rbh1_rbh2_fit)
      names(out) <- c("crbh.pairs", "crbh1", "crbh2", "rbh1_rbh2_fit")
      attr(out, "CRBHits.class") <- "crbh"
      attr(out, "crbh") <- TRUE
      attr(out, "keepSingleDirection") <- FALSE
      attr(out, "selfblast") <- selfblast
      return(out)
    }
    #if keepSingleDirection - include single - done
    if(keepSingleDirection){
      crbh1 <- data.frame(Map(c, cbind(rbh1, "rbh"), cbind(rbh1.sec[, 1:15], "sec"), cbind(single1[, 1:15], "single")))
      colnames(crbh1)[16] <- "rbh_class"
      crbh2 <- data.frame(Map(c, cbind(rbh2, "rbh"), cbind(rbh2.sec[, 1:15], "sec"), cbind(single2[, 1:15], "single")))
      colnames(crbh2)[16] <- "rbh_class"
      crbh <- data.frame(Map(c, crbh1[, c(1:2, 16)], crbh2[crbh2$rbh_class == "single", , drop = FALSE][, c(2, 1, 16)]))
      colnames(crbh) <- c("aa1", "aa2", "rbh_class")
      out <- list(crbh, crbh1, crbh2, rbh1_rbh2_fit)
      names(out) <- c("crbh.pairs", "crbh1", "crbh2", "rbh1_rbh2_fit")
      attr(out, "CRBHits.class") <- "crbh"
      attr(out, "crbh") <- TRUE
      attr(out, "keepSingleDirection") <- TRUE
      attr(out, "selfblast") <- selfblast
      return(out)
    }
  }
}
