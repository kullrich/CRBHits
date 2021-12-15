#' @title isoform2longest
#' @name isoform2longest
#' @description This function extracts the longest isoform from either NCBI or
#' ENSEMBL CDS input.
#' @param cds \code{DNAStringSet} [mandatory]
#' @param source source indicating either NCBI or ENSEMBL [default: NCBI]
#' @return \code{DNAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @importFrom curl curl_download
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' \dontrun{
#' ## load example sequence data
#' ## set NCBI URL
#' NCBI <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
#' HOMSAP.NCBI.cds.url <- paste0(NCBI,
#'     "GCF/000/001/405/GCF_000001405.39_GRCh38.p13/",
#'     "GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna.gz")
#' HOMSAP.NCBI.cds.file <- tempfile()
#' ## download CDS file
#' download.file(HOMSAP.NCBI.cds.url, HOMSAP.NCBI.cds.file, quiet=FALSE)
#' HOMSAP.NCBI.cds <- Biostrings::readDNAStringSet(HOMSAP.NCBI.cds.file)
#' ## get longest isoform
#' HOMSAP.NCBI.cds.longest <- isoform2longest(HOMSAP.NCBI.cds, "NCBI")
#' length(HOMSAP.NCBI.cds)
#' length(HOMSAP.NCBI.cds.longest)
#' ## set ENSEMBL URL
#' ensembl <- "ftp://ftp.ensembl.org/pub/release-101/fasta/"
#' ## set Homo sapiens CDS URL
#' HOMSAP.ENSEMBL.cds.url <- paste0(ensembl,
#'     "homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz")
#' HOMSAP.ENSEMBL.cds.file <- tempfile()
#' ## download CDS file
#' download.file(HOMSAP.ENSEMBL.cds.url, HOMSAP.ENSEMBL.cds.file, quiet=FALSE)
#' HOMSAP.ENSEMBL.cds <- Biostrings::readDNAStringSet(HOMSAP.ENSEMBL.cds.file)
#' ## get longest isoform
#' HOMSAP.ENSEMBL.cds.longest<-isoform2longest(HOMSAP.ENSEMBL.cds, "ENSEMBL")
#' length(HOMSAP.ENSEMBL.cds)
#' length(HOMSAP.ENSEMBL.cds.longest)
#' }
#' @export isoform2longest
#' @author Kristian K Ullrich

isoform2longest <- function(cds, source="NCBI"){
    if(source=="NCBI"){
        #get longer unique idx
        gene.idx.len <- length(grep("gene\\=", names(cds)))
        locus_tag.idx.len <- length(grep("locus_tag\\=", names(cds)))
        if(gene.idx.len >= locus_tag.idx.len){
            gene.id <- stringr::word(
            unlist(lapply(strsplit(names(cds)," \\[gene\\="),
                function(x) x[2])), 1)
        } else {
            gene.id <- stringr::word(
            unlist(lapply(strsplit(names(cds)," \\[locus_tag\\="),
                function(x) x[2])), 1)
        }
    }
    if(source=="ENSEMBL"){
        isoform.id <- stringr::word(names(cds), 4)
        isoform.id <- gsub("(\\D+)(\\.)(\\D{1})","\\1\\3", isoform.id)
        gene.id<-unlist(lapply(strsplit(isoform.id, "\\."), function(x) x[1]))
    }
    isoform.df <- data.frame(pos=seq(from=1, to=length(cds)),
    gene.width=Biostrings::width(cds),
          gene.id = gene.id)
    isoform.df.ordered <- isoform.df[order(isoform.df$gene.id,
        isoform.df$gene.width,
        decreasing = TRUE), ]
    longest.df <- isoform.df.ordered[!duplicated(isoform.df.ordered$gene.id), ]
    longest.idx <- sort(longest.df$pos)
    out <- cds[longest.idx]
    return(out)
}
