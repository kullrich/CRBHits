#' @title gtf2longest
#' @name gtf2longest
#' @description This function extracts the gene position from GTF input and
#' optional extracts the longest isoform.
#' @param gtffile \code{GTF path} [mandatory]
#' @param cds \code{DNAStringSet} [optional]
#' @param removeNonCoding specify if NonCoding transcripts should be removed
#' @param source source indicating either NCBI or ENSEMBL [default: NCBI]
#' @return \code{list}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word str_split_fixed str_split
#' @importFrom tidyr %>% unite
#' @importFrom dplyr filter select mutate summarise arrange distinct ungroup
#' desc
#' @importFrom readr read_tsv
#' @importFrom curl curl_download
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' \dontrun{
#' ## load example sequence data
#' ## set NCBI GTF URL
#' NCBI <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
#' ARATHA.NCBI.gtf.url <- paste0(NCBI,
#'     "GCF/000/001/735/GCF_000001735.4_TAIR10.1/",
#'     "GCF_000001735.4_TAIR10.1_genomic.gtf.gz")
#' ARATHA.NCBI.gtf.file <- tempfile()
#' ## download GTF file
#' download.file(ARATHA.NCBI.gtf.url, ARATHA.NCBI.gtf.file, quiet=FALSE)
#' ## set NCBI CDS URL
#' ARATHA.NCBI.cds.url <- paste0(NCBI,
#'     "GCF/000/001/735/GCF_000001735.4_TAIR10.1/",
#'     "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz")
#' ARATHA.NCBI.cds.file <- tempfile()
#' ## download CDS file
#' download.file(ARATHA.NCBI.cds.url, ARATHA.NCBI.cds.file, quiet=FALSE)
#' ## load CDS
#' ARATHA.NCBI.cds <- Biostrings::readDNAStringSet(ARATHA.NCBI.cds.file)
#' ## get genepos and longest isoform
#' ARATHA.NCBI.gtf.longest <- gtf2longest(gtffile=ARATHA.NCBI.gtf.file,
#'     cds=ARATHA.NCBI.cds, source="NCBI")
#' ARATHA.NCBI.gtf.longest$genepos
#' ARATHA.NCBI.gtf.longest$cds
#' ## set ENSEMBL GTF URL
#' ensembl <- "http://ftp.ensemblgenomes.org/pub/plants/release-52/"
#' ARATHA.ENSEMBL.gtf.url <- paste0(ensembl,
#'     "gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz")
#' ARATHA.ENSEMBL.gtf.file <- tempfile()
#' ## download GTF file
#' download.file(ARATHA.ENSEMBL.gtf.url, ARATHA.ENSEMBL.gtf.file, quiet=FALSE)
#' ## set ENSEMBL CDS URL
#' ARATHA.ENSEMBL.cds.url <- paste0(ensembl,
#'     "fasta/arabidopsis_thaliana/cds/",
#'     "Arabidopsis_thaliana.TAIR10.cds.all.fa.gz")
#' ARATHA.ENSEMBL.cds.file <- tempfile()
#' ## download CDS file
#' download.file(ARATHA.ENSEMBL.cds.url, ARATHA.ENSEMBL.cds.file, quiet=FALSE)
#' ARATHA.ENSEMBL.cds <- Biostrings::readDNAStringSet(ARATHA.ENSEMBL.cds.file)
#' ## get genepos and longest isoform
#' ARATHA.ENSEMBL.gtf.longest <- gtf2longest(gtffile=ARATHA.ENSEMBL.gtf.file,
#'     cds=ARATHA.ENSEMBL.cds)
#' ARATHA.ENSEMBL.gtf.longest$genepos
#' ARATHA.ENSEMBL.gtf.longest$cds
#' }
#' @export gtf2longest
#' @author Kristian K Ullrich

gtf2longest <- function(gtffile,
    cds=NULL,
    removeNonCoding=TRUE,
    source="NCBI"
    ){
    seqname <- NULL
    #source <- NULL
    feature <- NULL
    start <- NULL
    end <- NULL
    #score <- NULL
    strand <- NULL
    #frame <- NULL
    attribute <- NULL
    geneID <- NULL
    cdsID <- NULL
    transcriptID <- NULL
    transcriptLENGTH <- NULL
    mid <- NULL
    gtf <- readr::read_tsv(gtffile, col_names = FALSE, comment = "#")
    colnames(gtf) <- c("seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute")
    if(source=="NCBI"){
        ## extract gene
        gtf.gene <- gtf %>% dplyr::filter(feature == "gene")
        ## extract CDS
        gtf.CDS <- gtf %>% dplyr::filter(feature == "CDS")
        ## extract gene ID
        gtf.gene <- (gtf.gene %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("gene_id ", "",
            stringr::word(attribute, 1, sep=";"))))))
        gtf.CDS <- (gtf.CDS %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("gene_id ", "",
            stringr::word(attribute, 1, sep=";"))))))
        ## extract protein ID as transcriptID
        proteinID.idx <- grep(" protein_id", gtf.CDS$attribute)
        gtf.CDS <- gtf.CDS[proteinID.idx, ]
        proteinID <- stringr::str_split(gtf.CDS$attribute, ";")
        proteinID <- gsub(" ", "", gsub("\"", "",
            gsub("protein_id ", "",
            unlist(lapply(proteinID, function(x) x[grep(" protein_id", x)])))))
        gtf.CDS <- (gtf.CDS %>%
        dplyr::mutate(transcriptID = proteinID))
        ## extract width per isoform
        gtf.transcript <- (gtf.CDS %>%
            dplyr::group_by(transcriptID) %>%
                dplyr::summarise(
                transcriptID=unique(transcriptID),
                seqname=unique(seqname),
                start=min(start),
                end=max(end),
                strand=unique(strand),
                geneID=unique(geneID),
                transcriptLENGTH=sum(end-start+1)))
        ## retain only longest mRNA isoform
        gtf.transcript.longest <- (gtf.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            transcriptID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gtf.transcript.longest <- (gtf.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        ## add gene idx
        gtf.transcript.longest.ordered <- gtf.transcript.longest[
            order(gtf.transcript.longest$seqname,
            gtf.transcript.longest$mid,
            gtf.transcript.longest$transcriptID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gtf.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gtf.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gtf.transcript.longest.ordered <- dplyr::mutate(
            gtf.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gtf.transcript.genepos <- as.data.frame(
            gtf.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=transcriptID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx))
        attr(gtf.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            new_names <- unlist(lapply(stringr::str_split(
                unlist(lapply(
                stringr::str_split(stringr::word(names(cds)), "_cds_"),
                function(x) rev(x)[[1]])), "_"),
                function(x) paste0(rev(rev(x)[-1]), collapse="_")))
            new_names[which(new_names=="")]<-stringr::word(names(cds))[which(new_names=="")]
            names(cds) <- new_names
            cds <- cds[names(cds) %in% gtf.transcript.genepos$gene.seq.id]
            cds <- cds[match(gtf.transcript.genepos$gene.seq.id,
                stringr::word(names(cds)))]
        }
        return(setNames(list(gtf.transcript.genepos, cds),
            c("genepos", "cds")))
    }
    if(source=="ENSEMBL"){
        ## extract gene
        gtf.gene <- gtf %>% dplyr::filter(feature == "gene")
        ## extract transcript
        gtf.transcript <- gtf %>% dplyr::filter(feature == "transcript")
        ## extract CDS
        gtf.CDS <- gtf %>% dplyr::filter(feature == "CDS")
        ## extract gene ID
        gtf.gene <- (gtf.gene %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("gene_id ", "",
            stringr::word(attribute, 1, sep=";"))))))
        gtf.transcript <- (gtf.transcript %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("gene_id ", "",
            stringr::word(attribute, 1, sep=";"))))))
        gtf.CDS <- (gtf.CDS %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("gene_id ", "",
            stringr::word(attribute, 1, sep=";"))))))
        ## extract transcript ID and transcript VERSION
        gtf.transcript.attribute <- stringr::str_split(
            gtf.transcript$attribute, ";", simplify=TRUE)
        gtf.transcriptID <- apply(gtf.transcript.attribute, 1,
            function(x){
                gsub(" ", "", gsub("\"", "",
                gsub("transcript_id ", "" ,
                x[grep("transcript_id", x)])))
            })
        gtf.transcriptVERSION <- apply(gtf.transcript.attribute, 1,
            function(x){
                gsub(" ", "", gsub("\"", "",
                gsub("transcript_version ", "" ,
                x[grep("transcript_version", x)])))
            })
        if(length(gtf.transcriptVERSION)!=0){
            gtf.transcript <- (gtf.transcript %>%
            dplyr::mutate(transcriptID = paste0(gtf.transcriptID,
                ".", gtf.transcriptVERSION)))
        }
        if(length(gtf.transcriptVERSION)==0){
            gtf.transcript <- (gtf.transcript %>%
            dplyr::mutate(transcriptID = gtf.transcriptID))
        }
        gtf.CDS.attribute <- stringr::str_split(
            gtf.CDS$attribute, ";", simplify=TRUE)
        gtf.CDS.transcriptID <- apply(gtf.CDS.attribute, 1,
            function(x){
                gsub(" ", "", gsub("\"", "",
                gsub("transcript_id ", "" ,
                x[grep("transcript_id", x)])))
            })
        gtf.CDS.transcriptVERSION <- apply(gtf.CDS.attribute, 1,
            function(x){
                gsub(" ", "", gsub("\"", "",
                gsub("transcript_version ", "" ,
                x[grep("transcript_version", x)])))
            })
        if(length(gtf.CDS.transcriptVERSION)!=0){
            gtf.CDS <- (gtf.CDS %>%
            dplyr::mutate(transcriptID = paste0(gtf.CDS.transcriptID,
                ".", gtf.CDS.transcriptVERSION)))
        }
        if(length(gtf.CDS.transcriptVERSION)==0){
            gtf.CDS <- (gtf.CDS %>%
            dplyr::mutate(transcriptID = gtf.CDS.transcriptID))
        }
        ## extract width per isoform
        gtf.CDS.len <- (gtf.CDS %>%
        dplyr::group_by(transcriptID) %>%
            dplyr::summarise(transcriptID=unique(transcriptID),
            len=sum(end-start+1)))
        gtf.transcript <- (gtf.transcript %>%
        dplyr::mutate(transcriptLENGTH=
            gtf.CDS.len$len[
                match(gtf.transcript$transcriptID,
                gtf.CDS.len$transcriptID)]))
        ## retain only longest mRNA isoform
        gtf.transcript.longest <- (gtf.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            transcriptID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gtf.transcript.longest <- (gtf.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        ## add gene idx
        gtf.transcript.longest.ordered <- gtf.transcript.longest[
            order(gtf.transcript.longest$seqname,
            gtf.transcript.longest$mid,
            gtf.transcript.longest$transcriptID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gtf.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gtf.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gtf.transcript.longest.ordered <- dplyr::mutate(
            gtf.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gtf.transcript.genepos <- as.data.frame(
            gtf.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=transcriptID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx))
        attr(gtf.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            cds <- cds[stringr::word(names(cds)) %in%
                gtf.transcript.genepos$gene.seq.id]
            cds <- cds[match(gtf.transcript.genepos$gene.seq.id
                [gtf.transcript.genepos$gene.seq.id %in%
                    stringr::word(names(cds))],
                stringr::word(names(cds)))]
        }
        return(setNames(list(gtf.transcript.genepos, cds),
            c("genepos", "cds")))
    }
}
