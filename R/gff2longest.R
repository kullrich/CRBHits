#' @title gff2longest
#' @name gff2longest
#' @description This function extracts the gene position from GFF3 input and
#' optional extracts the longest isoform.
#' @param gff3file \code{GFF3 path} [mandatory]
#' @param cds \code{DNAStringSet} [optional]
#' @param removeNonCoding specify if NonCoding transcripts should be removed
#' @param source source indicating either NCBI or ENSEMBL [default: NCBI]
#' @return \code{list}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word str_split_fixed str_split
#' @importFrom tidyr %>% unite
#' @importFrom dplyr filter select mutate summarise arrange distinct ungroup
#' desc reframe
#' @importFrom readr read_tsv
#' @importFrom curl curl_download
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' \dontrun{
#' ## load example sequence data
#' ## set NCBI GFF3 URL
#' NCBI <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
#' ARATHA.NCBI.gff3.url <- paste0(NCBI,
#'     "GCF/000/001/735/GCF_000001735.4_TAIR10.1/",
#'     "GCF_000001735.4_TAIR10.1_genomic.gff.gz")
#' ARATHA.NCBI.gff3.file <- tempfile()
#' ## download GTF file
#' download.file(ARATHA.NCBI.gff3.url, ARATHA.NCBI.gff3.file, quiet=FALSE)
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
#' ARATHA.NCBI.gff3.longest <- gff2longest(gtffile=ARATHA.NCBI.gff3.file,
#'     cds=ARATHA.NCBI.cds, source="NCBI")
#' ARATHA.NCBI.gff3.longest$genepos
#' ARATHA.NCBI.gff3.longest$cds
#' ## set ENSEMBL GFF3 URL
#' ensembl <- "http://ftp.ensemblgenomes.org/pub/plants/release-52/"
#' ARATHA.ENSEMBL.gff3.url <- paste0(ensembl,
#'     "gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gff3.gz")
#' ARATHA.ENSEMBL.gff3.file <- tempfile()
#' ## download GTF file
#' download.file(ARATHA.ENSEMBL.gff3.url, ARATHA.ENSEMBL.gff3.file, quiet=FALSE)
#' ## set ENSEMBL CDS URL
#' ARATHA.ENSEMBL.cds.url <- paste0(ensembl,
#'     "fasta/arabidopsis_thaliana/cds/",
#'     "Arabidopsis_thaliana.TAIR10.cds.all.fa.gz")
#' ARATHA.ENSEMBL.cds.file <- tempfile()
#' ## download CDS file
#' download.file(ARATHA.ENSEMBL.cds.url, ARATHA.ENSEMBL.cds.file, quiet=FALSE)
#' ARATHA.ENSEMBL.cds <- Biostrings::readDNAStringSet(ARATHA.ENSEMBL.cds.file)
#' ## get genepos and longest isoform
#' ARATHA.ENSEMBL.gff3.longest <- gff2longest(gff3file=ARATHA.ENSEMBL.gff3.file,
#'     cds=ARATHA.ENSEMBL.cds)
#' ARATHA.ENSEMBL.gff3.longest$genepos
#' ARATHA.ENSEMBL.gff3.longest$cds
#' }
#' @export gff2longest
#' @author Kristian K Ullrich

gff2longest <- function(gff3file,
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
    if(source=="NCBI"){
        gff3 <- readr::read_tsv(gff3file, col_names = FALSE, comment = "#")
        colnames(gff3) <- c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute")
        rm_idx <- c(which(gff3$end<gff3$start), which(gff3$strand=="?"))
        if(length(rm_idx)>0){gff3 <- gff3[-rm_idx, ]}
        ## extract gene
        gff3.gene <- gff3 %>% dplyr::filter(feature == "gene")
        ## extract mRNA
        gff3.mRNA <- gff3 %>% dplyr::filter(feature == "mRNA")
        ## extract CDS
        gff3.CDS <- gff3 %>% dplyr::filter(feature == "CDS")
        ## extract ID
        gff3.gene <- (gff3.gene %>%
        dplyr::mutate(dbxref_geneID = unlist(lapply(lapply(
            stringr::str_split(gff3.gene$attribute, ";"), function(x) {
                    x[grep("Dbxref=", x)]}), function(x) {
                    stringr::str_split_fixed(x, "GeneID:", 2)[,2]}))))
        gff3.gene <- (gff3.gene %>%
        dplyr::mutate(geneID = unlist(lapply(lapply(
            stringr::str_split(gff3.gene$attribute, ";"), function(x) {
                    x[grep("ID=", x)]}), function(x) {
                    stringr::str_split_fixed(x, "ID=", 2)[,2]}))))
        gff3.gene <- (gff3.gene %>%
        dplyr::mutate(locusTAG = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.gene$attribute, ";"), function(x) {
            x[grep("locus_tag=", x)]}), function(x) {
            stringr::str_split_fixed(x, "locus_tag=", 2)[,2]}),
            function(x) {
            ifelse(length(x)==1, x, "")}))))
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(dbxref_geneID = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.mRNA$attribute, ";"), function(x) {
            x[grep("Dbxref=", x)]}), function(x) {
            unlist(stringr::str_split(x, ","))}), function(x) {
            stringr::str_split_fixed(x[grep("GeneID:", x)],
            "GeneID:", 2)[,2]}))))
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(transcriptID = unlist(lapply(lapply(
            stringr::str_split(gff3.mRNA$attribute, ";"), function(x) {
            x[grep("ID=", x)]}), function(x) {
            stringr::str_split_fixed(x, "ID=", 2)[,2]}))))
        ## extract mRNA$parentID == gene$geneID
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(geneID = unlist(lapply(lapply(
            stringr::str_split(gff3.mRNA$attribute, ";"), function(x) {
            x[grep("Parent=", x)]}), function(x) {
            stringr::str_split_fixed(x, "Parent=", 2)[,2]}))))
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(locusTAG = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.mRNA$attribute, ";"), function(x) {
            x[grep("locus_tag=", x)]}), function(x) {
            stringr::str_split_fixed(x, "locus_tag=", 2)[,2]}),
            function(x) {
            ifelse(length(x)==1, x, "")}))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(dbxref_geneID = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.CDS$attribute, ";"), function(x) {
            x[grep("Dbxref=", x)]}), function(x) {
            unlist(stringr::str_split(x, ","))}), function(x) {
            stringr::str_split_fixed(x[grep("GeneID:", x)],
            "GeneID:", 2)[,2]}))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(cdsID = unlist(lapply(lapply(
            stringr::str_split(gff3.CDS$attribute, ";"), function(x) {
            x[grep("ID=", x)]}), function(x) {
            stringr::str_split_fixed(x, "ID=", 2)[,2]}))))
        ## extract CDS$parentID == mRNA$transcriptID
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(transcriptID = unlist(lapply(lapply(
            stringr::str_split(gff3.CDS$attribute, ";"), function(x) {
            x[grep("Parent=", x)]}), function(x) {
            stringr::str_split_fixed(x, "Parent=", 2)[,2]}))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(proteinID = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.CDS$attribute, ";"), function(x) {
            x[grep("protein_id=", x)]}), function(x) {
            stringr::str_split_fixed(x, "protein_id=", 2)[,2]}),
            function(x) {
            ifelse(length(x)==1, x, "")}))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(locusTAG = unlist(lapply(lapply(lapply(
            stringr::str_split(gff3.CDS$attribute, ";"), function(x) {
            x[grep("locus_tag=", x)]}), function(x) {
            stringr::str_split_fixed(x, "locus_tag=", 2)[,2]}),
            function(x) {
            ifelse(length(x)==1, x, "")}))))
        ## add geneID
        gff3.CDS <- (gff3.CDS %>%
            dplyr::mutate(geneID =
            gff3.mRNA$geneID[match(gff3.CDS$transcriptID,
            gff3.mRNA$transcriptID)]))
        ## if geneID is NA use transcriptID
        gff3.CDS[which(is.na(gff3.CDS$geneID)),]$geneID <- gff3.CDS[which(
            is.na(gff3.CDS$geneID)),]$transcriptID
        ## extract width per isoform
        gff3.transcript <- (gff3.CDS %>%
            dplyr::group_by(transcriptID) %>%
                dplyr::reframe(
                transcriptID=unique(transcriptID),
                seqname=unique(seqname),
                start=min(start),
                end=max(end),
                strand=unique(strand),
                geneID=unique(geneID),
                cdsID=unique(cdsID),
                proteinID=unique(proteinID),
                dbxref_geneID=unique(dbxref_geneID),
                transcriptLENGTH=sum(end-start+1)))
        ## retain only longest mRNA isoform
        gff3.transcript.longest <- (gff3.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            cdsID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gff3.transcript.longest <- (gff3.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        ## add gene idx
        gff3.transcript.longest.ordered <- gff3.transcript.longest[
            order(gff3.transcript.longest$seqname,
            gff3.transcript.longest$mid,
            gff3.transcript.longest$cdsID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gff3.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gff3.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gff3.transcript.longest.ordered <- dplyr::mutate(
            gff3.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gff3.transcript.genepos <- as.data.frame(
            gff3.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=cdsID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx,
                geneID=geneID,
                dbxref_geneID=dbxref_geneID,
                transcriptID=transcriptID,
                cdsID=cdsID,
                proteinID=proteinID))
        gff3.transcript.genepos$gene.seq.id <- gsub("cds-", "",
            gff3.transcript.genepos$gene.seq.id)
        attr(gff3.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            new_names <- unlist(lapply(stringr::str_split(
                unlist(lapply(
                stringr::str_split(stringr::word(names(cds)), "_cds_"),
                function(x) rev(x)[[1]])), "_"),
                function(x) paste0(rev(rev(x)[-1]), collapse="_")))
            new_names[which(new_names=="")]<-stringr::word(names(cds))[which(new_names=="")]
            names(cds) <- new_names
            cds <- cds[names(cds) %in% gff3.transcript.genepos$gene.seq.id]
            cds <- cds[match(gff3.transcript.genepos$gene.seq.id,
                stringr::word(names(cds)))]
        }
        return(setNames(list(gff3.transcript.genepos, cds),
            c("genepos", "cds")))
    }
    if(source=="ENSEMBL"){
        gff3 <- readLines(gff3file)
        gff3 <- gff3[substr(gff3,1,1)!="#"]
        gff3 <- as.data.frame(stringr::str_split_fixed(gff3, "\t", 9))
        colnames(gff3) <- c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute")
        gff3$start <- as.integer(gff3$start)
        gff3$end <- as.integer(gff3$end)
        ## extract gene
        gff3.gene <- gff3 %>% dplyr::filter(feature == "gene")
        ## extract mRNA
        gff3.mRNA <- gff3 %>% dplyr::filter(feature == "mRNA")
        ## extract CDS
        gff3.CDS <- gff3 %>% dplyr::filter(feature == "CDS")
        ## extract ID
        gff3.gene <- (gff3.gene %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(transcriptID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(cdsID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        ## extract geneID
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("Parent=", "",
            stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        ## extract transcriptID
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(transcriptID = gsub(" ", "", gsub("\"", "",
            gsub("Parent=", "",
            stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        ## add geneID
        gff3.CDS <- (gff3.CDS %>%
            dplyr::mutate(geneID =
            gff3.mRNA$geneID[match(gff3.CDS$transcriptID,
            gff3.mRNA$transcriptID)]))
        ## if geneID is NA use transcriptID
        gff3.CDS[which(is.na(gff3.CDS$geneID)),]$geneID <- gff3.CDS[which(
            is.na(gff3.CDS$geneID)),]$transcriptID
        ## extract width per isoform
        gff3.transcript <- (gff3.CDS %>%
            dplyr::group_by(transcriptID) %>%
                dplyr::summarise(
                transcriptID=unique(transcriptID),
                seqname=unique(seqname),
                start=min(start),
                end=max(end),
                strand=unique(strand),
                geneID=unique(geneID),
                cdsID=unique(cdsID),
                transcriptLENGTH=sum(end-start+1)))
        ## retain only longest mRNA isoform
        gff3.transcript.longest <- (gff3.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            cdsID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gff3.transcript.longest <- (gff3.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        ## add gene idx
        gff3.transcript.longest.ordered <- gff3.transcript.longest[
            order(gff3.transcript.longest$seqname,
            gff3.transcript.longest$mid,
            gff3.transcript.longest$cdsID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gff3.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gff3.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gff3.transcript.longest.ordered <- dplyr::mutate(
            gff3.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gff3.transcript.genepos <- as.data.frame(
            gff3.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=cdsID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx))
        gff3.transcript.genepos$gene.seq.id <- gsub("CDS:", "",
            gff3.transcript.genepos$gene.seq.id)
        attr(gff3.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            cds <- cds[stringr::word(names(cds)) %in%
                gff3.transcript.genepos$gene.seq.id]
            cds <- cds[match(gff3.transcript.genepos$gene.seq.id,
                stringr::word(names(cds)))]
        }
        return(setNames(list(gff3.transcript.genepos, cds),
            c("genepos", "cds")))
    }
    if(source=="GEMOMA"){
        gff3 <- read.table(gff3file, sep="\t", header=FALSE)
        gff3 <- as.data.frame(gff3)
        colnames(gff3) <- c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute")
        gff3$start <- as.integer(gff3$start)
        gff3$end <- as.integer(gff3$end)
        gff3.gene <- gff3 %>% dplyr::filter(feature == "gene")
        gff3.mRNA <- gff3 %>% dplyr::filter(feature == "mRNA")
        gff3.CDS <- gff3 %>% dplyr::filter(feature == "CDS")
        ## extract ID
        gff3.gene <- (gff3.gene %>%
        dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(transcriptID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        ## extract geneID
        gff3.mRNA <- (gff3.mRNA %>%
        dplyr::mutate(geneID = stringr::str_split_fixed(gsub("=", "",
            stringr::str_split_fixed(gff3.mRNA$attribute,
                "Parent", 2)[,2]), ";", 2)[,1]))
        ## extract transcriptID
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(transcriptID = stringr::str_split_fixed(gsub("=", "",
            stringr::str_split_fixed(gff3.CDS$attribute,
                "Parent", 2)[,2]), ";", 2)[,1]))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(geneID = gff3.mRNA$geneID[match(gff3.CDS$transcriptID,gff3.mRNA$transcriptID)]))
        ## if geneID is NA use transcriptID
        gff3.CDS[which(is.na(gff3.CDS$geneID)),]$geneID <- gff3.CDS[which(
            is.na(gff3.CDS$geneID)),]$transcriptID
        ## extract width per isoform
        gff3.transcript <- (gff3.CDS %>%
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
        gff3.transcript.longest <- (gff3.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            transcriptID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gff3.transcript.longest <- (gff3.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        gff3.transcript.longest.ordered <- gff3.transcript.longest[
            order(gff3.transcript.longest$seqname,
            gff3.transcript.longest$mid,
            gff3.transcript.longest$transcriptID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gff3.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gff3.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gff3.transcript.longest.ordered <- dplyr::mutate(
            gff3.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gff3.transcript.genepos <- as.data.frame(
            gff3.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=transcriptID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx))
        attr(gff3.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            cds <- cds[stringr::word(names(cds)) %in%
                gff3.transcript.genepos$gene.seq.id]
        }
        return(setNames(list(gff3.transcript.genepos, cds),
            c("genepos", "cds")))
    }
    if(source=="BAKTA"){
        gff3 <- readr::read_tsv(gff3file, col_names = FALSE, comment = "#")
        colnames(gff3) <- c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute")
        ## extract gene
        #gff3.gene <- gff3 %>% dplyr::filter(feature == "gene")
        ## extract mRNA
        #gff3.mRNA <- gff3 %>% dplyr::filter(feature == "mRNA")
        ## extract CDS
        gff3.CDS <- gff3 %>% dplyr::filter(feature == "CDS")
        ## extract ID
        #gff3.gene <- (gff3.gene %>%
        #dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
        #    gsub("ID=", "",
        #    stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        #gff3.mRNA <- (gff3.mRNA %>%
        #dplyr::mutate(transcriptID = gsub(" ", "", gsub("\"", "",
        #    gsub("ID=", "",
        #    stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(cdsID = gsub(" ", "", gsub("\"", "",
            gsub("ID=", "",
            stringr::str_split_fixed(attribute, ";", 2)[,1])))))
        ## extract geneID
        #gff3.mRNA <- (gff3.mRNA %>%
        #dplyr::mutate(geneID = gsub(" ", "", gsub("\"", "",
        #    gsub("Parent=", "",
        #    stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        ## extract transcriptID
        #gff3.CDS <- (gff3.CDS %>%
        #dplyr::mutate(transcriptID = gsub(" ", "", gsub("\"", "",
        #    gsub("Parent=", "",
        #    stringr::str_split_fixed(attribute, ";", 3)[,2])))))
        ## add geneID
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(transcriptID = gff3.CDS$cdsID))
        #gff3.CDS <- (gff3.CDS %>%
        #    dplyr::mutate(geneID =
        #    gff3.mRNA$geneID[match(gff3.CDS$transcriptID,
        #    gff3.mRNA$transcriptID)]))
        gff3.CDS <- (gff3.CDS %>%
        dplyr::mutate(geneID = gff3.CDS$cdsID))
        ## if geneID is NA use transcriptID
        #gff3.CDS[which(is.na(gff3.CDS$geneID)),]$geneID <- gff3.CDS[which(
        #    is.na(gff3.CDS$geneID)),]$transcriptID
        ## extract width per isoform
        gff3.transcript <- (gff3.CDS %>%
            dplyr::group_by(transcriptID) %>%
                dplyr::summarise(
                transcriptID=unique(transcriptID),
                seqname=unique(seqname),
                start=min(start),
                end=max(end),
                strand=unique(strand),
                geneID=unique(geneID),
                cdsID=unique(cdsID),
                transcriptLENGTH=sum(end-start+1)))
        ## retain only longest mRNA isoform
        gff3.transcript.longest <- (gff3.transcript %>%
            dplyr::arrange(seqname, desc(transcriptLENGTH),
            cdsID, geneID, start) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(geneID, .keep_all=TRUE) %>%
            dplyr::mutate(mid=start + (end-start)/2))
        if(removeNonCoding){
            gff3.transcript.longest <- (gff3.transcript.longest %>%
                dplyr::filter(!is.na(transcriptLENGTH)))
        }
        ## add gene idx
        gff3.transcript.longest.ordered <- gff3.transcript.longest[
            order(gff3.transcript.longest$seqname,
            gff3.transcript.longest$mid,
            gff3.transcript.longest$cdsID,
            decreasing=FALSE, method="radix"), , drop=FALSE]
        gene.mid.idx <- gff3.transcript.longest.ordered %>%
            dplyr::select(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- gff3.transcript.longest.ordered %>%
            dplyr::distinct(seqname, mid) %>%
            tidyr::unite("gene.chr.mid", sep=":")
        gene.mid.idx.red <- dplyr::mutate(gene.mid.idx.red,
            gene.idx=seq(from=1, to=dim(gene.mid.idx.red)[1]))
        gene.idx <- gene.mid.idx.red$gene.idx[match(gene.mid.idx$gene.chr.mid,
            gene.mid.idx.red$gene.chr.mid)]
        gff3.transcript.longest.ordered <- dplyr::mutate(
            gff3.transcript.longest.ordered, gene.idx=gene.idx)
        ## create gene position matrix to be used for downstream analysis
        gff3.transcript.genepos <- as.data.frame(
            gff3.transcript.longest.ordered %>%
            dplyr::select(
                gene.seq.id=cdsID,
                gene.chr=seqname,
                gene.start=start,
                gene.end=end,
                gene.mid=mid,
                gene.strand=strand,
                gene.idx=gene.idx))
        gff3.transcript.genepos$gene.seq.id <- gsub("cds-", "",
            gff3.transcript.genepos$gene.seq.id)
        attr(gff3.transcript.genepos, "CRBHits.class") <- "genepos"
        if(!is.null(cds)){
            names(cds) <- apply(stringr::str_split_fixed(
                unlist(lapply(
                stringr::str_split(stringr::word(names(cds)), "_cds_"),
                function(x) rev(x)[[1]])), "_", 3)[, c(1, 2)], 1,
                function(x) paste0(x, collapse="_"))
            cds <- cds[names(cds) %in% gff3.transcript.genepos$gene.seq.id]
            cds <- cds[match(gff3.transcript.genepos$gene.seq.id,
                stringr::word(names(cds)))]
        }
        return(setNames(list(gff3.transcript.genepos, cds),
            c("genepos", "cds")))
    }
}
