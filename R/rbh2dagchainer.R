#' @title rbh2dagchainer
#' @name rbh2dagchainer
#' @description This function runs DAGchainer (http://dagchainer.sourceforge.net/) given CRBHit pairs and gene positions for both cds1 and cds2. The default options are set to not compare gene positions in base pairs but instead using gene order (gene.idx).
#' @param crbh see(\code{\link[CRBHits]{cds2rbh}} [mandatory]
#' @param gene.position.cds1 specify gene position for cds1 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param gene.position.cds2 specify gene position for cds2 sequences (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param plotDotPlot specify if dotplot should be plotted [default: FALSE]
#' @param dagchainerpath specify the PATH to the DAGchainer binaries [default: /extdata/dagchainer/]
#' @param gap_open_penalty gap open penalty [default: 0]
#' @param gap_extension_penalty gap extension penalty [default: -3]
#' @param gap_length length of a gap (avgerage distance expected between two syntenic genes); if type is set to "bp" use 10000 [default: 1]
#' @param max_match_score Maximum match score [default: 50]
#' @param max_dist_allowed maximum distance allowed between two matches; if type is set to "bp" use 200000 [default: 200000]
#' @param max_evalue Maximum E-value [default: 1e-3]
#' @param min_number_aligned_pairs Minimum number of Aligned Pairs [default: 5]
#' @param type specify if gene order index "idx" or gene base pair position "bp" should be extracted and used with DAGchainer [default: idx]
#' @return \code{DAGchanier} results if type == "bp"\cr
#' 1: $gene1.chr\cr
#' 2: $gene1.seq.id\cr
#' 3: $gene1.start\cr
#' 4: $gene1.end\cr
#' 5: $gene1.mid\cr
#' 6: $gene2.chr\cr
#' 7: $gene2.seq.id\cr
#' 8: $gene2.start\cr
#' 9: $gene2.end\cr
#' 10: $gene2.mid\cr
#' 11: $evalue\cr
#' 12: $score\cr
#' \cr
#' \code{DAGchanier} results if type == "idx"\cr
#' 1: $gene1.chr\cr
#' 2: $gene1.seq.id\cr
#' 3: $gene1.idx1\cr
#' 4: $gene1.idx2\cr
#' 5: $gene1.idx\cr
#' 6: $gene2.chr\cr
#' 7: $gene2.seq.id\cr
#' 8: $gene2.idx1\cr
#' 9: $gene2.idx2\cr
#' 10: $gene2.idx\cr
#' 11: $evalue\cr
#' 12: $score\cr
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @export rbh2dagchainer
#' @author Kristian K Ullrich

rbh2dagchainer <- function(crbh,
                           gene.position.cds1 = NULL,
                           gene.position.cds2 = NULL,
                           plotDotPlot = FALSE,
                           dagchainerpath = paste0(find.package("CRBHits"),
                                             "/extdata/dagchainer/"),
                           gap_open_penalty = 0,
                           gap_extension_penalty = -3,
                           gap_length = 1,
                           max_match_score = 50,
                           max_dist_allowed = 10,
                           max_evalue = 1e-3,
                           min_number_aligned_pairs = 5,
                           type = "idx"
                           ){
  if(!dir.exists(dagchainerpath)){stop("Error: DAGchainer PATH does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.dagchainer() function.")}
  if(!file.exists(paste0(dagchainerpath, "dagchainer"))){stop("Error: dagchainer binary does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.dagchainer() function.")}
  if(!file.exists(paste0(dagchainerpath, "run_DAG_chainer.pl"))){stop("Error: run_DAG_chainer.pl does not exist. Please specify correct PATH and/or look into package installation prerequisites. Try to use make.dagchainer() function.")}
  genepos.colnames <- c("gene.seq.id", "gene.chr", "gene.start", "gene.end",
                        "gene.mid", "gene.strand", "gene.idx")
  if(!is.null(gene.position.cds1) & !is.null(gene.position.cds2)){
    genepos.colnames <- c("gene.seq.id", "gene.chr", "gene.start", "gene.end",
                          "gene.mid", "gene.strand", "gene.idx")
    if(any(colnames(gene.position.cds1) != genepos.colnames) | any(colnames(gene.position.cds2) != genepos.colnames)){
      stop("Error: Please specify gene.position.cds1 and gene.position.cds2 as indicated in cds2genepos()")
    }
  }
  if(type == "bp"){
    aa1.genepos <- gene.position.cds1[match(crbh$crbh.pairs$aa1,
                                            gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
      dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
    aa2.genepos <- gene.position.cds2[match(crbh$crbh.pairs$aa2,
                                            gene.position.cds2$gene.seq.id), , drop =FALSE] %>%
      dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
    dagchainer.input <- cbind(aa1.genepos, aa2.genepos)
    colnames(dagchainer.input) <- c("gene1.chr", "gene1.seq.id", "gene1.start", "gene1.end",
                                    "gene2.chr", "gene2.seq.id", "gene2.start", "gene2.end")
    dagchainer.evalue <- crbh$crbh1[match(dagchainer.input$gene1.seq.id, crbh$crbh1$query_id), , drop =FALSE]$evalue
    dagchainer.input <- cbind(dagchainer.input, evalue = dagchainer.evalue)
  }
  if(type == "idx"){
    aa1.genepos <- gene.position.cds1[match(crbh$crbh.pairs$aa1,
                                            gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
      dplyr::select(gene.chr, gene.seq.id, gene.idx1 = gene.idx, gene.idx2 = gene.idx)
    aa2.genepos <- gene.position.cds2[match(crbh$crbh.pairs$aa2,
                                            gene.position.cds2$gene.seq.id), , drop =FALSE] %>%
      dplyr::select(gene.chr, gene.seq.id, gene.idx1 = gene.idx, gene.idx2 = gene.idx)
    dagchainer.input <- cbind(aa1.genepos, aa2.genepos)
    colnames(dagchainer.input) <- c("gene1.chr", "gene1.seq.id", "gene1.idx1", "gene1.idx2",
                                    "gene2.chr", "gene2.seq.id", "gene2.idx1", "gene2.idx2")
    dagchainer.evalue <- crbh$crbh1[match(dagchainer.input$gene1.seq.id, crbh$crbh1$query_id), , drop =FALSE]$evalue
    dagchainer.input <- cbind(dagchainer.input, evalue = dagchainer.evalue)
  }
  tmp <- tempfile()
  write.table(dagchainer.input,
              sep="\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE,
              file = tmp)
  system(paste0(dagchainerpath, "run_DAG_chainer.pl -i ", tmp,
                " -o ", sprintf("%0i",gap_open_penalty),
                " -e ", sprintf("%0i",gap_extension_penalty),
                " -g ", sprintf("%0i",gap_length),
                " -M ", sprintf("%0i",max_match_score),
                " -D ", sprintf("%0i",max_dist_allowed),
                " -E ", max_evalue,
                " -A ", sprintf("%0i",min_number_aligned_pairs)), ignore.stdout = TRUE, ignore.stderr = TRUE)
  dagchainer.results <- read.table(paste0(tmp, ".aligncoords"), sep = "\t", header = FALSE)
  if(type == "bp"){
    colnames(dagchainer.results) <- c("gene1.chr", "gene1.seq.id", "gene1.start", "gene1.end",
                                      "gene2.chr", "gene2.seq.id", "gene2.start", "gene2.end",
                                      "evalue", "score")
    #add gene1.mid and gene2.mid for plotting
    dagchainer.results <- dagchainer.results %>%
      dplyr::mutate(gene1.mid = (gene1.start + gene1.end)/2,
                    gene2.mid = (gene2.start + gene2.end)/2) %>%
      dplyr::select(gene1.chr, gene1.seq.id, gene1.start, gene1.end, gene1.mid,
                    gene2.chr, gene2.seq.id, gene2.start, gene2.end, gene2.mid,
                    evalue, score)
  }
  if(type == "idx"){
    colnames(dagchainer.results) <- c("gene1.chr", "gene1.seq.id", "gene1.idx1", "gene1.idx2",
                                      "gene2.chr", "gene2.seq.id", "gene2.idx1", "gene2.idx2",
                                      "evalue", "score")
    #add gene1.mid and gene2.mid for plotting
    dagchainer.results <- dagchainer.results %>%
      dplyr::mutate(gene1.idx = gene1.idx1,
                    gene2.idx = gene2.idx1) %>%
      dplyr::select(gene1.chr, gene1.seq.id, gene1.idx1, gene1.idx2, gene1.idx,
                    gene2.chr, gene2.seq.id, gene2.idx1, gene2.idx2, gene2.idx,
                    evalue, score)
  }
  if(plotDotPlot){

  }
  return(dagchainer.results)
}
