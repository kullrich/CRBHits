#' @title rbh2dagchainer
#' @name rbh2dagchainer
#' @description This function runs DAGchainer
#' (http://dagchainer.sourceforge.net/) given CRBHit pairs and gene positions
#' for both cds1 and cds2. The default options are set to not compare gene
#' positions in base pairs but instead using gene order (gene.idx).
#' @param rbhpairs (conditional-)reciprocal best hit (CRBHit) pair result
#' (see \code{\link[CRBHits]{cds2rbh}}) [mandatory]
#' @param selfblast1 (conditional-)reciprocal best hit (CRBHit) pair selfblast
#' result for cds1 sequences (see \code{\link[CRBHits]{cds2rbh}}) [optional]
#' @param selfblast2 (conditional-)reciprocal best hit (CRBHit) pair selfblast
#' result for cds2 sequences (see \code{\link[CRBHits]{cds2rbh}}) [optional]
#' @param gene.position.cds1 specify gene position for cds1 sequences
#' (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param gene.position.cds2 specify gene position for cds2 sequences
#' (see \code{\link[CRBHits]{cds2genepos}}) [default: NULL]
#' @param dagchainerpath specify the PATH to the DAGchainer binaries
#' [default: /extdata/dagchainer/]
#' @param gap_open_penalty gap open penalty [default: 0]
#' @param gap_extension_penalty gap extension penalty [default: -3]
#' @param gap_length length of a gap (avgerage distance expected between two
#' syntenic genes); if type is set to "idx" use 1 [default: 10000]
#' @param max_match_score Maximum match score [default: 50]
#' @param max_dist_allowed maximum distance allowed between two matches;
#' if type is set to "idx" use 20 [default: 200000]
#' @param max_evalue Maximum E-value [default: 1e-3]
#' @param ignore_tandem ignore tandem duplicates [default = TRUE]
#' @param only_tandem only tandem alignments [default = FALSE]
#' @param min_number_aligned_pairs Minimum number of Aligned Pairs [default: 5]
#' @param type specify if gene order index "idx" or gene base pair position
#' "bp" should be extracted and used with DAGchainer [default: bp]
#' @param plotDotPlot specify if dotplot should be plotted [default: FALSE]
#' @param DotPlotTitle specify DotPlot title [default: 'DAGchainer results']
#' @param colorBy specify if dagchainer groups should be colored by "Ka", "Ks",
#' "Ka/Ks" or "none" [default: none]
#' @param kaks specify Ka/Ks input obtained via `rbh2kaks()` [default: NULL]
#' @param ka.max specify max Ka to be filtered [default: 5]
#' @param ks.max specify max Ks to be filtered [default: 5]
#' @param ka.min specify min Ka to be filtered [default: 0]
#' @param ks.min specify min Ks to be filtered [default: 0]
#' @param select.chr filter results for chromosome names [default: NULL]
#' @return \code{DAGchanier} results\cr
#' 1: $gene1.chr\cr
#' 2: $gene1.seq.id\cr
#' 3: $gene1.start\cr
#' 4: $gene1.end\cr
#' 5: $gene1.mid\cr
#' 6: $gene1.idx\cr
#' 7: $gene2.chr\cr
#' 8: $gene2.seq.id\cr
#' 9: $gene2.start\cr
#' 10: $gene2.end\cr
#' 11: $gene2.mid\cr
#' 12: $gene2.idx\cr
#' 13: $evalue\cr
#' 14: $score\cr
#' @importFrom tidyr %>%
#' @importFrom dplyr bind_cols select group_by group_map group_keys mutate
#' @importFrom stringr word
#' @importFrom ggplot2 ggplot geom_point geom_abline facet_wrap
#' @importFrom utils write.table
#' @seealso \code{\link[CRBHits]{plot_dagchainer}},
#' \code{\link[CRBHits]{cds2genepos}},
#' \code{\link[CRBHits]{tandemdups}},
#' \code{\link[CRBHits]{rbh2kaks}}
#' @references Haas BJ et al. (2004) DAGchainer: a tool for mining segmental
#' genome duplications and synteny. \emph{Bioinformatics.}
#' \bold{20(18)}, 3643-3646.
#' @examples
#' ## compile dagchainer
#' CRBHits::make_dagchainer()
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' ## get selfhits CRBHit pairs
#' ath_selfhits_crbh <- cds2rbh(
#'     cds1=ath,
#'     cds2=ath,
#'     plotCurve=TRUE)
#' ## get gene position
#' ath.genepos <- cds2genepos(
#'     cds=ath,
#'     source="ENSEMBL")
#' ## get DAGchainer results
#' ath_selfblast_crbh.dagchainer <- rbh2dagchainer(
#'     rbhpairs=ath_selfhits_crbh,
#'     gene.position.cds1=ath.genepos,
#'     gene.position.cds2=ath.genepos)
#' head(ath_selfblast_crbh.dagchainer)
#' ## plot dagchainer
#' plot_dagchainer(
#'     dag=ath_selfblast_crbh.dagchainer)
#' @export rbh2dagchainer
#' @author Kristian K Ullrich

rbh2dagchainer <- function(rbhpairs,
    selfblast1=NULL,
    selfblast2=NULL,
    gene.position.cds1=NULL,
    gene.position.cds2=NULL,
    dagchainerpath=paste0(find.package("CRBHits"),
        "/extdata/dagchainer/"),
    gap_open_penalty=0,
    gap_extension_penalty=-3,
    gap_length=10000,
    max_match_score=50,
    max_dist_allowed=200000,
    max_evalue=1e-3,
    ignore_tandem=TRUE,
    only_tandem=FALSE,
    min_number_aligned_pairs=5,
    type="bp",
    plotDotPlot=FALSE,
    DotPlotTitle="DAGchainer results",
    colorBy="none",
    kaks=NULL,
    ka.max=5,
    ks.max=5,
    ka.min=0,
    ks.min=0,
    select.chr=NULL){
    gene.seq.id <- NULL
    gene1.seq.id <- NULL
    gene2.seq.id <- NULL
    gene.chr <- NULL
    gene1.chr <- NULL
    gene2.chr <- NULL
    gene.start <- NULL
    gene.end <- NULL
    gene.idx <- NULL
    gene1.idx1 <- NULL
    gene1.idx2 <- NULL
    gene2.idx1 <- NULL
    gene2.idx2 <- NULL
    evalue <- NULL
    score <- NULL
    if(ignore_tandem == TRUE & only_tandem == TRUE){
        stop("both options set to TRUE can not go together. Set one of them to
            FALSE")
    }
    if(attributes(rbhpairs)$CRBHits.class != "crbh"){
        stop("Please obtain rbhpairs via the cds2rbh or the cdsfile2rbh
            function")
    }
    if(is.null(gene.position.cds1) | is.null(gene.position.cds2)){
        stop("Please provide gene position information. Can be obtained via the
            cds2genepos() function")
    }
    if(!is.null(gene.position.cds1)){
        if(attributes(gene.position.cds1)$CRBHits.class != "genepos"){
            stop("Please obtain gene position via the cds2genepos() function or
                add a 'genepos' class attribute")
        }
    }
    if(!is.null(gene.position.cds2)){
        if(attributes(gene.position.cds2)$CRBHits.class != "genepos"){
            stop("Please obtain gene position via the cds2genepos() function or
                add a 'genepos' class attribute")
        }
    }
    if(!dir.exists(dagchainerpath)){
        stop("Error: DAGchainer PATH does not exist. Please specify correct
            PATH and/or look into package installation prerequisites. Try to
            use make_dagchainer() function.")
    }
    if(!file.exists(paste0(dagchainerpath, "dagchainer"))){
        stop("Error: dagchainer binary does not exist. Please specify correct
            PATH and/or look into package installation prerequisites. Try to
            use make_dagchainer() function.")
    }
    if(!file.exists(paste0(dagchainerpath, "run_DAG_chainer.pl"))){
        stop("Error: run_DAG_chainer.pl does not exist. Please specify correct
            PATH and/or look into package installation prerequisites. Try to use
            make_dagchainer() function.")
    }
    selfblast <- attributes(rbhpairs)$selfblast
    genepos.colnames <- c("gene.seq.id", "gene.chr", "gene.start", "gene.end",
        "gene.mid", "gene.strand", "gene.idx")
    if(plotDotPlot){
        if(colorBy == "Ka" | colorBy == "Ks" | colorBy == "Ka/Ks"){
            if(is.null(kaks)){
                stop("Please provide Ka/Ks input as obtained via rbh2kaks()")
            }
            if(!is.null(kaks)){
                if(attributes(kaks)$CRBHits.class != "kaks"){
                    stop("Please obtain Ka/Ks input via the rbh2kaks() function
                        or add a 'kaks' class attribute")
                }
            }
        }
    }
    if(type=="bp"){
        aa1.genepos <- gene.position.cds1[match(rbhpairs$crbh.pairs$aa1,
            gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
            dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
        aa2.genepos <- gene.position.cds2[match(rbhpairs$crbh.pairs$aa2,
            gene.position.cds2$gene.seq.id), , drop =FALSE] %>%
            dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
        dagchainer.input <- cbind(aa1.genepos, aa2.genepos)
        colnames(dagchainer.input) <- c("gene1.chr", "gene1.seq.id",
            "gene1.start", "gene1.end", "gene2.chr", "gene2.seq.id",
            "gene2.start", "gene2.end")
        dagchainer.evalue <- rbhpairs$crbh1[match(dagchainer.input$gene1.seq.id,
            rbhpairs$crbh1$query_id), , drop=FALSE]$evalue
        dagchainer.input <- cbind(dagchainer.input, evalue=dagchainer.evalue)
        if(selfblast){
            dagchainer.input$gene1.chr<-paste0("AA1:",
                dagchainer.input$gene1.chr)
            dagchainer.input$gene2.chr<-paste0("AA1:",
                dagchainer.input$gene2.chr)
        } else{
            dagchainer.input$gene1.chr<-paste0("AA1:",
                dagchainer.input$gene1.chr)
            dagchainer.input$gene2.chr<-paste0("AA2:",
                dagchainer.input$gene2.chr)
        }
        if(!is.null(selfblast1)){
            selfblast1.aa1.genepos <- gene.position.cds1[
                match(selfblast1$crbh.pairs$aa1,
                gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
            selfblast1.aa2.genepos <- gene.position.cds1[
                match(selfblast1$crbh.pairs$aa2,
                gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
            selfblast1.dagchainer.input <- cbind(selfblast1.aa1.genepos,
                selfblast1.aa2.genepos)
            colnames(selfblast1.dagchainer.input) <- c("gene1.chr",
                "gene1.seq.id", "gene1.start", "gene1.end", "gene2.chr",
                "gene2.seq.id", "gene2.start", "gene2.end")
            selfblast1.dagchainer.evalue <- selfblast1$crbh1[
                match(selfblast1.dagchainer.input$gene1.seq.id,
                selfblast1$crbh1$query_id), , drop=FALSE]$evalue
            selfblast1.dagchainer.input <- cbind(selfblast1.dagchainer.input,
                evalue=selfblast1.dagchainer.evalue)
            selfblast1.dagchainer.input$gene1.chr<-paste0("AA1:",
                selfblast1.dagchainer.input$gene1.chr)
            selfblast1.dagchainer.input$gene2.chr<-paste0("AA1:",
                selfblast1.dagchainer.input$gene2.chr)
            dagchainer.input <- rbind(dagchainer.input,
                selfblast1.dagchainer.input)
        }
        if(!is.null(selfblast2)){
            selfblast2.aa1.genepos <- gene.position.cds2[
                match(selfblast2$crbh.pairs$aa1,
                gene.position.cds2$gene.seq.id), , drop=FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
            selfblast2.aa2.genepos <- gene.position.cds2[
                match(selfblast2$crbh.pairs$aa2,
                gene.position.cds2$gene.seq.id), , drop=FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.start, gene.end)
            selfblast2.dagchainer.input <- cbind(selfblast2.aa1.genepos,
                    selfblast2.aa2.genepos)
            colnames(selfblast2.dagchainer.input) <- c("gene1.chr",
                "gene1.seq.id", "gene1.start", "gene1.end", "gene2.chr",
                "gene2.seq.id", "gene2.start", "gene2.end")
            selfblast2.dagchainer.evalue <- selfblast2$crbh1[
                match(selfblast2.dagchainer.input$gene1.seq.id,
                selfblast2$crbh1$query_id), , drop=FALSE]$evalue
            selfblast2.dagchainer.input <- cbind(selfblast2.dagchainer.input,
                evalue=selfblast2.dagchainer.evalue)
            selfblast2.dagchainer.input$gene1.chr<-paste0("AA2:",
                selfblast2.dagchainer.input$gene1.chr)
            selfblast2.dagchainer.input$gene2.chr<-paste0("AA2:",
                selfblast2.dagchainer.input$gene2.chr)
            dagchainer.input <- rbind(dagchainer.input,
                selfblast2.dagchainer.input)
        }
    }
    if(type=="idx"){
        aa1.genepos <- gene.position.cds1[match(rbhpairs$crbh.pairs$aa1,
            gene.position.cds1$gene.seq.id), , drop=FALSE] %>%
        dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
        gene.idx2=gene.idx)
        aa2.genepos <- gene.position.cds2[match(rbhpairs$crbh.pairs$aa2,
            gene.position.cds2$gene.seq.id), , drop=FALSE] %>%
        dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
            gene.idx2=gene.idx)
        dagchainer.input <- cbind(aa1.genepos, aa2.genepos)
        colnames(dagchainer.input) <- c("gene1.chr", "gene1.seq.id",
            "gene1.idx1", "gene1.idx2", "gene2.chr", "gene2.seq.id",
            "gene2.idx1", "gene2.idx2")
        dagchainer.evalue <- rbhpairs$crbh1[match(dagchainer.input$gene1.seq.id,
            rbhpairs$crbh1$query_id), , drop=FALSE]$evalue
        dagchainer.input <- cbind(dagchainer.input, evalue=dagchainer.evalue)
        if(selfblast){
            dagchainer.input$gene1.chr<-paste0("AA1:",
                dagchainer.input$gene1.chr)
            dagchainer.input$gene2.chr<-paste0("AA1:",
                dagchainer.input$gene2.chr)
        } else{
            dagchainer.input$gene1.chr<-paste0("AA1:",
                dagchainer.input$gene1.chr)
            dagchainer.input$gene2.chr<-paste0("AA2:",
                dagchainer.input$gene2.chr)
        }
        if(!is.null(selfblast1)){
            selfblast1.aa1.genepos <- gene.position.cds1[
                match(selfblast1$crbh.pairs$aa1,
                gene.position.cds1$gene.seq.id), , drop =FALSE] %>%
            dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
            gene.idx2=gene.idx)
            selfblast1.aa2.genepos <- gene.position.cds1[
                match(selfblast1$crbh.pairs$aa2,
                gene.position.cds1$gene.seq.id), , drop=FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
                gene.idx2=gene.idx)
            selfblast1.dagchainer.input <- cbind(selfblast1.aa1.genepos,
                selfblast1.aa2.genepos)
            colnames(selfblast1.dagchainer.input) <- c("gene1.chr",
                "gene1.seq.id", "gene1.idx1", "gene1.idx2", "gene2.chr",
                "gene2.seq.id", "gene2.idx1", "gene2.idx2")
            selfblast1.dagchainer.evalue <- selfblast1$crbh1[
                match(selfblast1.dagchainer.input$gene1.seq.id,
                selfblast1$crbh1$query_id), , drop=FALSE]$evalue
            selfblast1.dagchainer.input <- cbind(selfblast1.dagchainer.input,
                evalue=selfblast1.dagchainer.evalue)
            selfblast1.dagchainer.input$gene1.chr<-paste0("AA1:",
                selfblast1.dagchainer.input$gene1.chr)
            selfblast1.dagchainer.input$gene2.chr<-paste0("AA1:",
                selfblast1.dagchainer.input$gene2.chr)
            dagchainer.input <- rbind(dagchainer.input,
                selfblast1.dagchainer.input)
        }
        if(!is.null(selfblast2)){
            selfblast2.aa1.genepos <- gene.position.cds2[
                match(selfblast2$crbh.pairs$aa1,
                gene.position.cds2$gene.seq.id), , drop=FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
                gene.idx2=gene.idx)
            selfblast2.aa2.genepos <- gene.position.cds2[
                match(selfblast2$crbh.pairs$aa2,
                gene.position.cds2$gene.seq.id), , drop=FALSE] %>%
                dplyr::select(gene.chr, gene.seq.id, gene.idx1=gene.idx,
                gene.idx2=gene.idx)
            selfblast2.dagchainer.input <- cbind(selfblast2.aa1.genepos,
                selfblast2.aa2.genepos)
            colnames(selfblast2.dagchainer.input) <- c("gene1.chr",
                "gene1.seq.id", "gene1.idx1", "gene1.idx2", "gene2.chr",
                "gene2.seq.id", "gene2.idx1", "gene2.idx2")
            selfblast2.dagchainer.evalue <- selfblast2$crbh1[
                match(selfblast2.dagchainer.input$gene1.seq.id,
                selfblast2$crbh1$query_id), , drop=FALSE]$evalue
            selfblast2.dagchainer.input <- cbind(selfblast2.dagchainer.input,
                evalue=selfblast2.dagchainer.evalue)
            selfblast2.dagchainer.input$gene1.chr<-paste0("AA2:",
                selfblast2.dagchainer.input$gene1.chr)
            selfblast2.dagchainer.input$gene2.chr<-paste0("AA2:",
                selfblast2.dagchainer.input$gene2.chr)
            dagchainer.input <- rbind(dagchainer.input,
                selfblast2.dagchainer.input)
        }
    }
    tmp <- tempfile()
    write.table(dagchainer.input,
        sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=tmp)
    dagchainercmd <- paste0(dagchainerpath, "run_DAG_chainer.pl")
    dagchainerargs <- c(
        "-i", tmp,
        "-o", sprintf("%0i", gap_open_penalty),
        "-e", sprintf("%0i", gap_extension_penalty),
        "-g", sprintf("%0i", gap_length),
        "-M", sprintf("%0i", max_match_score),
        "-D", sprintf("%0i", max_dist_allowed),
        "-E", max_evalue,
        "-A", sprintf("%0i", min_number_aligned_pairs))
    if(selfblast | (!is.null(selfblast1) | !is.null(selfblast2))){
        if(ignore_tandem && !only_tandem){
            print(paste(dagchainercmd,
                paste(dagchainerargs, collapse = " "), "-I", "-s"))
            system2(command=dagchainercmd, args=c(dagchainerargs, "-I", "-s"),
                stdout=FALSE, stderr=FALSE)
        }
        if(only_tandem && !ignore_tandem){
            print(paste(dagchainercmd,
                paste(dagchainerargs, collapse = " "), "-T", "-s"))
            system2(command=dagchainercmd, args=c(dagchainerargs, "-T", "-s"),
                stdout=FALSE, stderr=FALSE)
        }
        if(!ignore_tandem && !only_tandem){
            print(paste(dagchainercmd,
                paste(dagchainerargs, collapse = " "), "-s"))
            system2(command=dagchainercmd, args=c(dagchainerargs, "-s"),
                stdout=FALSE, stderr=FALSE)
        }
    } else {
        print(paste(dagchainercmd,
            paste(dagchainerargs, collapse = " ")))
        system2(command=dagchainercmd, args=dagchainerargs,
            stdout=FALSE, stderr=FALSE)
    }
    if(file.info(paste0(tmp, ".aligncoords"))[["size"]]==0){
        return(NULL)
    }
    dagchainer.results <- read.table(paste0(tmp, ".aligncoords"), sep = "\t",
        header = FALSE)
    dagchainer.groups <- readLines(paste0(tmp, ".aligncoords"))
    dagchainer.groups <- dagchainer.groups[grep("num aligned pairs",
        dagchainer.groups)]
    dagchainer.groups.len <- as.numeric(stringr::str_split_fixed(
        stringr::str_split_fixed(dagchainer.groups,
        "num aligned pairs: ", 2)[, 2],"\\)", 2)[, 1])
    dagchainer.groups.chr1 <- stringr::str_split_fixed(
        stringr::str_split_fixed(dagchainer.groups,
        "alignment ", 2)[, 2]," ", 2)[, 1]
    dagchainer.groups.chr2 <- stringr::str_split_fixed(
        stringr::str_split_fixed(dagchainer.groups,
        "vs. ", 2)[, 2]," ", 2)[, 1]
    dagchainer.groups.idx <- stringr::str_split_fixed(
        stringr::str_split_fixed(dagchainer.groups,
        "Alignment #", 2)[, 2]," ", 2)[, 1]
    dagchainer.groups.idx.reverse <- grep("(reverse)", dagchainer.groups)
    dagchainer.groups.id <- paste0(dagchainer.groups.chr1,
        ":", dagchainer.groups.chr2, ":", dagchainer.groups.idx)
    dagchainer.groups.id[dagchainer.groups.idx.reverse] <- paste0(
        dagchainer.groups.id[dagchainer.groups.idx.reverse], ".rev")
    dagchainer.groups.out <- unlist(apply(cbind(dagchainer.groups.id,
        dagchainer.groups.len), 1, function(x) rep(x[1], x[2])))
    #if(!selfblast){
    #  dagchainer.groups.out <- gsub("AA1:", "", dagchainer.groups.out)
    #  dagchainer.groups.out <- gsub("AA2:", "", dagchainer.groups.out)
    #}
    if(type=="bp"){
        colnames(dagchainer.results) <- c("gene1.chr", "gene1.seq.id",
            "gene1.start", "gene1.end", "gene2.chr", "gene2.seq.id",
            "gene2.start", "gene2.end", "evalue", "score")
        #if(!selfblast){
        #  dagchainer.results$gene1.chr <- gsub("AA1:", "",
        #dagchainer.results$gene1.chr)
        #  dagchainer.results$gene2.chr <- gsub("AA2:", "",
        #dagchainer.results$gene2.chr)
        #}
        #add gene1.mid and gene2.mid for plotting
        dagchainer.results <- dagchainer.results %>%
            dplyr::mutate(gene1.mid=gene1.start+((gene1.end - gene1.start)/2),
            gene2.mid=gene2.start+((gene2.end - gene2.start)/2)) %>%
            dplyr::select(gene1.chr, gene1.seq.id, gene1.start, gene1.end,
            gene1.mid, gene2.chr, gene2.seq.id, gene2.start, gene2.end,
            gene2.mid, evalue, score)
        gene1.idx <- rep(NA, length(dagchainer.results$gene1.chr))
        gene2.idx <- rep(NA, length(dagchainer.results$gene2.chr))
        #process AA1
        gene1.chr.AA1.pos <- which(substr(
            dagchainer.results$gene1.chr, 1, 3)=="AA1")
        gene2.chr.AA1.pos <- which(substr(
            dagchainer.results$gene2.chr, 1, 3)=="AA1")
        if(length(gene1.chr.AA1.pos)>0){
            gene1.idx[gene1.chr.AA1.pos] <- gene.position.cds1$gene.idx[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
        }
        if(length(gene2.chr.AA1.pos)>0){
            gene2.idx[gene2.chr.AA1.pos] <- gene.position.cds1$gene.idx[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
        }
        #process AA2
        gene1.chr.AA2.pos <- which(substr(
            dagchainer.results$gene1.chr, 1, 3)=="AA2")
        gene2.chr.AA2.pos <- which(substr(
            dagchainer.results$gene2.chr, 1, 3)=="AA2")
        if(length(gene1.chr.AA2.pos)>0){
            gene1.idx[gene1.chr.AA2.pos] <- gene.position.cds2$gene.idx[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
        }
        if(length(gene2.chr.AA2.pos)>0){
            gene2.idx[gene2.chr.AA2.pos] <- gene.position.cds2$gene.idx[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
        }
        dagchainer.results <- dagchainer.results %>%
            dplyr::mutate(gene1.idx = gene1.idx, gene2.idx = gene2.idx) %>%
            dplyr::select(gene1.chr, gene1.seq.id, gene1.start, gene1.end,
            gene1.mid, gene1.idx, gene2.chr, gene2.seq.id, gene2.start,
            gene2.end, gene2.mid, gene2.idx, evalue, score)
    }
    if(type=="idx"){
        colnames(dagchainer.results) <- c("gene1.chr", "gene1.seq.id",
            "gene1.idx1", "gene1.idx2", "gene2.chr", "gene2.seq.id",
            "gene2.idx1", "gene2.idx2", "evalue", "score")
        #if(!selfblast){
        #  dagchainer.results$gene1.chr <- gsub("AA1:", "",
        #dagchainer.results$gene1.chr)
        #  dagchainer.results$gene2.chr <- gsub("AA2:", "",
        #dagchainer.results$gene2.chr)
        #}
        #add gene1.mid and gene2.mid for plotting
        dagchainer.results <- dagchainer.results %>%
            dplyr::mutate(gene1.idx = gene1.idx1, gene2.idx = gene2.idx1) %>%
            dplyr::select(gene1.chr, gene1.seq.id, gene1.idx1, gene1.idx2,
            gene1.idx, gene2.chr, gene2.seq.id, gene2.idx1, gene2.idx2,
            gene2.idx, evalue, score)
            gene1.start <- rep(NA, length(dagchainer.results$gene1.chr))
            gene2.start <- rep(NA, length(dagchainer.results$gene2.chr))
            gene1.end <- rep(NA, length(dagchainer.results$gene1.chr))
            gene2.end <- rep(NA, length(dagchainer.results$gene2.chr))
            gene1.mid <- rep(NA, length(dagchainer.results$gene1.chr))
            gene2.mid <- rep(NA, length(dagchainer.results$gene2.chr))
            #process AA1
        gene1.chr.AA1.pos <- which(substr(
        dagchainer.results$gene1.chr, 1, 3)=="AA1")
        gene2.chr.AA1.pos <- which(substr(
        dagchainer.results$gene2.chr, 1, 3)=="AA1")
        if(length(gene1.chr.AA1.pos)>0){
            gene1.start[gene1.chr.AA1.pos] <- gene.position.cds1$gene.start[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
            gene1.end[gene1.chr.AA1.pos] <- gene.position.cds1$gene.end[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
            gene1.mid[gene1.chr.AA1.pos] <- gene.position.cds1$gene.mid[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
        }
        if(length(gene2.chr.AA1.pos)>0){
            gene2.start[gene2.chr.AA1.pos] <- gene.position.cds1$gene.start[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
            gene2.end[gene2.chr.AA1.pos] <- gene.position.cds1$gene.end[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
            gene2.mid[gene2.chr.AA1.pos] <- gene.position.cds1$gene.mid[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA1.pos],
                gene.position.cds1$gene.seq.id)]
        }
        #process AA2
        gene1.chr.AA2.pos <- which(substr(
            dagchainer.results$gene1.chr, 1, 3)=="AA2")
        gene2.chr.AA2.pos <- which(substr(
            dagchainer.results$gene2.chr, 1, 3)=="AA2")
        if(length(gene1.chr.AA2.pos)>0){
            gene1.start[gene1.chr.AA2.pos] <- gene.position.cds2$gene.start[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
            gene1.end[gene1.chr.AA2.pos] <- gene.position.cds2$gene.end[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
            gene1.mid[gene1.chr.AA2.pos] <- gene.position.cds2$gene.mid[
                match(dagchainer.results$gene1.seq.id[gene1.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
        }
        if(length(gene2.chr.AA2.pos)>0){
            gene2.start[gene2.chr.AA2.pos] <- gene.position.cds2$gene.start[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
            gene2.end[gene2.chr.AA2.pos] <- gene.position.cds2$gene.end[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
            gene2.mid[gene2.chr.AA2.pos] <- gene.position.cds2$gene.mid[
                match(dagchainer.results$gene2.seq.id[gene2.chr.AA2.pos],
                gene.position.cds2$gene.seq.id)]
        }
        dagchainer.results <- dagchainer.results %>%
            dplyr::mutate(gene1.start=gene1.start,
            gene2.start=gene2.start,
            gene1.end=gene1.end,
            gene2.end=gene2.end,
            gene1.mid=gene1.mid,
            gene2.mid=gene2.mid) %>%
            dplyr::select(gene1.chr, gene1.seq.id, gene1.start, gene1.end,
            gene1.mid, gene1.idx, gene2.chr, gene2.seq.id, gene2.start,
            gene2.end, gene2.mid, gene2.idx, evalue, score)
    }
    #add dagchainer group
    dagchainer.results <- dagchainer.results %>%
        dplyr::mutate(dagchainer_group=dagchainer.groups.out)
    attr(dagchainer.results, "CRBHits.class") <- "dagchainer"
    attr(dagchainer.results, "selfblast") <- selfblast
    if(!is.null(selfblast1) | !is.null(selfblast2)){
        attr(dagchainer.results, "selfblast") <- TRUE
    }
    if(plotDotPlot){
        print(plot_dagchainer(dagchainer.results,
            DotPlotTitle=DotPlotTitle,
            colorBy=colorBy,
            kaks=kaks,
            ka.max=ka.max,
            ks.max=ks.max,
            ka.min=ka.min,
            ks.min=ks.min,
            select.chr=select.chr))
    }
    return(dagchainer.results)
}
