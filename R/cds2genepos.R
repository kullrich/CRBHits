#' @title cds2genepos
#' @name cds2genepos
#' @description This function extracts the gene position from either NCBI or
#' ENSEMBL CDS input.
#' @param cds \code{DNAStringSet} [mandatory]
#' @param source source indicating either NCBI or ENSEMBL [default: NCBI]
#' @param keep.names vector indicating gene ids to be kept before chromosomal
#' position assignment [default: NULL]
#' @return \code{matrix}
#' 1: $gene.seq.id\cr
#' 2: $gene.chr\cr
#' 3: $gene.start\cr
#' 4: $gene.end\cr
#' 5: $gene.mid\cr
#' 6: $gene.strand\cr
#' 7: $gene.idx\cr
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @importFrom dplyr select distinct mutate filter
#' @importFrom tidyr %>% unite
#' @importFrom curl curl_download
#' @seealso \code{\link[CRBHits]{isoform2longest}}
#' @examples
#' ## load example sequence data
#' data("ath", package="CRBHits")
#' ath.genepos <- cds2genepos(
#'     cds=ath,
#'     source="ENSEMBL")
#' \dontrun{
#' ## load example sequence data
#' ## set EnsemblPlants URL
#' ensemblPlants <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/"
#' ## set Arabidopsis thaliana CDS URL
#' ARATHA.cds.url <- paste0(ensemblPlants,
#'     "arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz")
#' ARATHA.cds.file <- tempfile()
#' ## download CDS
#' download.file(ARATHA.cds.url, ARATHA.cds.file, quiet=FALSE)
#' ARATHA.cds <- Biostrings::readDNAStringSet(ARATHA.cds.file)
#' ## get gene position
#' ARATHA.cds.genepos <- cds2genepos(ARATHA.cds, "ENSEMBL")
#' }
#' @export cds2genepos
#' @author Kristian K Ullrich

cds2genepos <- function(cds,
    source="NCBI",
    keep.names=NULL
    ){
    seq.id <- stringr::word(names(cds), 1)
    if(source == "NCBI"){
        gene.location <- unlist(lapply(strsplit(names(cds), " \\[location="),
            function(x) strsplit(x[2], "\\]")[[1]][1]))
        gene.start <- as.numeric(gsub("<", "",
            unlist(lapply(strsplit(gene.location, "\\.\\."),
            function(x) strsplit(rev(unlist(
                strsplit(x[1], "\\(")))[[1]],",")[[1]][1]))))
        gene.end <- as.numeric(gsub(">", "", unlist(lapply(
            strsplit(gene.location, "\\.\\."),
            function(x) gsub("\\)", "", rev(
                strsplit(rev(x)[[1]],",")[[1]][1]))))))
        gene.mid <- gene.start + ((gene.end - gene.start) / 2)
        gene.chr <- gsub("lcl\\|", "", unlist(lapply(strsplit(seq.id, "_cds_"),
            function(x) x[1])))
        gene.strand <- rep("1", length(cds))
        gene.strand[grep("complement", gene.location)] <- "-1"
    }
    if(source == "ENSEMBL"){
        gene.location <- stringr::word(names(cds), 3)
        gene.start <- as.numeric(unlist(lapply(strsplit(gene.location, "\\:"),
            function(x) x[4])))
        gene.end <- as.numeric(unlist(lapply(strsplit(gene.location, "\\:"),
            function(x) x[5])))
        gene.mid <- gene.start + ((gene.end - gene.start) / 2)
        gene.strand <- unlist(lapply(strsplit(gene.location, "\\:"),
            function(x) x[6]))
        gene.chr <- unlist(lapply(strsplit(gene.location, "\\:"),
            function(x) x[3]))
    }
    gene.pos.df <- data.frame(gene.seq.id=seq.id,
        gene.chr=gene.chr,
        gene.start=gene.start,
        gene.end=gene.end,
        gene.mid=gene.mid,
        gene.strand=gene.strand)
    if(!is.null(keep.names)){
        gene.pos.df <- gene.pos.df[gene.pos.df$gene.seq.id %in% keep.names , ,
            drop=FALSE]
    }
    gene.pos.df.ordered <- gene.pos.df[order(gene.pos.df$gene.chr,
        gene.pos.df$gene.mid,
        gene.pos.df$gene.seq.id,
        decreasing=FALSE, method="radix"), , drop=FALSE]
    gene.mid.idx <- gene.pos.df.ordered %>%
        dplyr::select(gene.chr, gene.mid) %>%
        tidyr::unite("gene.chr.mid", sep=":")
    gene.mid.idx.red <- gene.pos.df.ordered %>%
        dplyr::distinct(gene.chr, gene.mid) %>%
        tidyr::unite("gene.chr.mid", sep=":")
    gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
        gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
    gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
        gene.mid.idx.red$gene.chr.mid)]
    gene.pos.df.ordered <- dplyr::mutate(gene.pos.df.ordered, gene.idx=gene.idx)
    attr(gene.pos.df.ordered, "CRBHits.class") <- "genepos"
    return(gene.pos.df.ordered)
}
