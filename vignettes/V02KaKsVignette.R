## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## load vignette specific libraries
library(CRBHits)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(curl))
## compile LAST, KaKs_Calculator2.0 and DAGchainer for the vignette
vignette.paths <- make_vignette()

## -----------------------------------------------------------------------------
## set URLs for Arabidopis thaliana and Arabidopsis lyrata from NCBI Genomes

## set NCBI URL
NCBI <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
## set Arabidopsis thaliana CDS URL
ARATHA.cds.url <- paste0(NCBI, "GCF/000/001/735/GCF_000001735.4_TAIR10.1/",
                      "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz")

## set Arabidopsis lyrata CDS URL
ARALYR.cds.url <- paste0(NCBI, "GCF/000/004/255/GCF_000004255.2_v.1.0/",
                      "GCF_000004255.2_v.1.0_cds_from_genomic.fna.gz")
## get Arabidopsis thaliana CDS
ARATHA.cds <- Biostrings::readDNAStringSet(ARATHA.cds.url)
## get Arabidopsis lyrata CDS
ARALYR.cds <- Biostrings::readDNAStringSet(ARALYR.cds.url)

## -----------------------------------------------------------------------------
## get longest isoforms

## get longest isoform
ARATHA.cds.longest <- isoform2longest(ARATHA.cds, "NCBI")
ARALYR.cds.longest <- isoform2longest(ARALYR.cds, "NCBI")

## get help
#?cds2longest

## -----------------------------------------------------------------------------
## calculate CRBHit pairs for A. thaliana and A. lyrata using 2 threads

## input from CDS obtained from NCBI
## longest isoform selection
## query coverage >= 50%
## rost199 filter
ARATHA_ARALYR_crbh <- cds2rbh(ARATHA.cds,
                              ARALYR.cds,
                              qcov = 0.5, rost1999 = TRUE,
                              longest.isoform = TRUE,
                              isoform.source = "NCBI",
                              threads = 2, plotCurve = TRUE,
                              lastpath = vignette.paths[1])

## get help
#?cds2rbh

## -----------------------------------------------------------------------------
## get gene position idx from NCBI CDS

## extract gene position from CDS
ARATHA.cds.genepos <- cds2genepos(ARATHA.cds, source = "NCBI")
## extract gene position from longest isoform CDS
ARATHA.cds.longest.genepos <- cds2genepos(ARATHA.cds.longest, source = "NCBI")
## show first entries
head(ARATHA.cds.genepos)
head(ARATHA.cds.longest.genepos)
## get number of gene isoforms with same index
table(table(ARATHA.cds.genepos$gene.idx))
table(table(ARATHA.cds.longest.genepos$gene.idx))

##get help
#?cds2genepos

## -----------------------------------------------------------------------------
## get gene position idx from GTF/GFF3

## set NCBI URL
ensemblPlants <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-48/"
## set Arabidopsis thaliana GFF3 URL
ARATHA.GFF.url <- paste0(ensemblPlants, "gff3/arabidopsis_thaliana/",
                      "Arabidopsis_thaliana.TAIR10.48.gff3.gz")
## downlaod and gunzip file
ARATHA.GFF.file <- tempfile()
ARATHA.GFF.file.gz <- paste0(ARATHA.GFF.file, ".gz")
download.file(ARATHA.GFF.url, ARATHA.GFF.file.gz, quiet = FALSE)
system(paste0("gunzip -f ", ARATHA.GFF.file.gz))
## read ARATHA.GFF.file
ARATHA.gff <- read.table(ARATHA.GFF.file, sep = "\t", 
                         quote = "", header = FALSE)
colnames(ARATHA.gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

## extract gene
ARATHA.gff.gene <- ARATHA.gff %>% dplyr::filter(feature == "gene")
## get gene ID
ARATHA.gff.gene.id <- gsub("ID\\=", "", 
                           stringr::word(ARATHA.gff.gene$attribute,
                           1, sep=";"))
## extract mRNA
ARATHA.gff.mRNA <- ARATHA.gff %>% dplyr::filter(feature == "mRNA")
## get mRNA ID
ARATHA.gff.mRNA.id <- gsub("ID\\=", "", 
                           stringr::word(ARATHA.gff.mRNA$attribute,
                           1, sep=";"))
## get mRNA based gene ID
ARATHA.gff.mRNA.parent.id <- gsub("Parent\\=", "", 
                                  stringr::word(ARATHA.gff.mRNA$attribute,
                                  2, sep=";"))
## add mRNA ID and mRNA based gene ID
ARATHA.gff.mRNA <- ARATHA.gff.mRNA %>% 
                     dplyr::mutate(gene.id = ARATHA.gff.mRNA.parent.id,
                                   mRNA.id = ARATHA.gff.mRNA.id)
## extract CDS
ARATHA.gff.CDS <- ARATHA.gff %>% dplyr::filter(feature == "CDS")
## get CDS based mRNA ID
ARATHA.gff.CDS.parent.id <- gsub("Parent\\=", "", 
                                 stringr::word(ARATHA.gff.CDS$attribute,
                                 2, sep=";"))
## add CDS based mRNA ID
ARATHA.gff.CDS <- ARATHA.gff.CDS %>%
                    dplyr::mutate(mRNA.id = ARATHA.gff.CDS.parent.id)
## get width per mRNA isoform
ARATHA.gff.mRNA.len <- ARATHA.gff.CDS %>% dplyr::group_by(mRNA.id) %>%
                         dplyr::summarise(mRNA.id = unique(mRNA.id),
                                          len = sum(end-start))
## add mRNA isoform width
ARATHA.gff.mRNA <- ARATHA.gff.mRNA %>% 
                    dplyr::mutate(mRNA.len = 
                    ARATHA.gff.mRNA.len$len[
                        match(ARATHA.gff.mRNA$mRNA.id, 
                              ARATHA.gff.mRNA.len$mRNA.id)])
## retain only longest mRNA isoform
ARATHA.gff.mRNA.longest <- ARATHA.gff.mRNA %>% dplyr::arrange(seqname, start,
  gene.id, desc(mRNA.len)) %>% dplyr::distinct(gene.id, .keep_all = TRUE) %>% 
  dplyr::mutate(gene.mid = (start+end)/2)
## add gene idx
ARATHA.gff.mRNA.longest <- ARATHA.gff.mRNA.longest %>% 
  dplyr::mutate(gene.idx = seq(from = 1, to = dim(ARATHA.gff.mRNA.longest)[1]))
## create gene position matrix to be used for downstream analysis
ARATHA.gff.genepos <- ARATHA.gff.mRNA.longest %>% dplyr::select(
  gene.seq.id = mRNA.id,
  gene.chr = seqname,
  gene.start = start,
  gene.end = end,
  gene.mid = gene.mid,
  gene.strand = strand,
  gene.idx = gene.idx)
## set attribute
attr(ARATHA.gff.genepos, "CRBHits.class") <- "genepos"
## show first entries of generated gene position
head(ARATHA.gff.genepos)

## -----------------------------------------------------------------------------
## example to assign tandem duplicates given selfblast CRBHit pairs and gene position

## get selfblast CRBHit pairs for A. thaliana
ARATHA_selfblast_crbh <- cds2rbh(ARATHA.cds,
                                 ARATHA.cds,
                                 qcov = 0.5, rost1999 = TRUE,
                                 longest.isoform = TRUE,
                                 isoform.source = "NCBI",
                                 threads = 2, plotCurve = TRUE,
                                 lastpath = vignette.paths[1])

## get selfblast CRBHit pairs for A. lyrata
ARALYR_selfblast_crbh <- cds2rbh(ARALYR.cds,
                                 ARALYR.cds,
                                 qcov = 0.5, rost1999 = TRUE,
                                 longest.isoform = TRUE,
                                 isoform.source = "NCBI",
                                 threads = 2, plotCurve = TRUE,
                                 lastpath = vignette.paths[1])

## get gene position for A. thaliana longest isoforms
ARATHA.cds.longest.genepos <- cds2genepos(ARATHA.cds.longest, source = "NCBI")
## get gene position for A. lyrata longest isoforms
ARALYR.cds.longest.genepos <- cds2genepos(ARALYR.cds.longest, source = "NCBI")

## assign tandem duplicates for A. thaliana
ARATHA.cds.longest.tandemdups <- tandemdups(rbhpairs = ARATHA_selfblast_crbh,
                                            genepos = ARATHA.cds.longest.genepos,
                                            dupdist = 5)

## assign tandem duplicates for A. lyrata
ARALYR.cds.longest.tandemdups <- tandemdups(rbhpairs = ARALYR_selfblast_crbh,
                                            genepos = ARALYR.cds.longest.genepos,
                                            dupdist = 5)
## get help
#?tandemdups

## -----------------------------------------------------------------------------
## example how to plot tandem duplicated gene groups

## get tandem group size
tandem_group_size <- ARATHA.cds.longest.tandemdups %>%
  dplyr::group_by(tandem_group) %>% dplyr::group_size()
table(tandem_group_size)

## use dplyr::mutate to assign group size
ARATHA.cds.longest.tandemdups <- ARATHA.cds.longest.tandemdups %>%
  dplyr::mutate(tandem_group_size = unlist(apply(cbind(tandem_group_size,
  tandem_group_size), 1, function(x) rep(x[1], x[2]))))

## use dplyr::group_by to plot group by chromosome and colored by group size
ARATHA.cds.longest.tandemdups %>% dplyr::group_by(gene.chr) %>%
  ggplot2::ggplot(aes(x = gene.mid, y = gene.mid)) +
  ggplot2::geom_point(shape = 20, aes(col = as.factor(tandem_group_size))) +
  ggplot2::facet_wrap(~ gene.chr) +
  ggplot2::scale_colour_manual(values = CRBHitsColors(length(table(tandem_group_size))))

## -----------------------------------------------------------------------------
## example how to run DAGchainer on CRBHit pairs and gene positions

## DAGchainer using gene base pair (gene bp start end)
## Note: change parameter to fit bp option
ARATHA_ARALYR_crbh.dagchainer.bp <- rbh2dagchainer(rbhpairs = ARATHA_ARALYR_crbh, 
                                      gene.position.cds1 = ARATHA.cds.longest.genepos,
                                      gene.position.cds2 = ARALYR.cds.longest.genepos,
                                      type = "bp",
                                      gap_length = 10000,
                                      max_dist_allowed = 200000,
                                      dagchainerpath = vignette.paths[3])

## DAGchainer using gene index (gene order)
ARATHA_ARALYR_crbh.dagchainer.idx <- rbh2dagchainer(rbhpairs = ARATHA_ARALYR_crbh, 
                                       gene.position.cds1 = ARATHA.cds.longest.genepos,
                                       gene.position.cds2 = ARALYR.cds.longest.genepos,
                                       type = "idx",
                                       gap_length = 1,
                                       max_dist_allowed = 20,
                                       dagchainerpath = vignette.paths[3])

## get help
#?rbh2dagchainer

## -----------------------------------------------------------------------------
## example how to plot pairwise chromosomal syntenic groups (DAGchainer results)

## plot DAGchainer results for each chromosome combination
plot_dagchainer(ARATHA_ARALYR_crbh.dagchainer.bp)

## plot DAGchainer results selected chromosomes
g <- plot_dagchainer(ARATHA_ARALYR_crbh.dagchainer.bp,
                     select.chr = c("NC_003070.9", "NC_003071.7", "NC_003074.8",
                                    "NC_003075.7", "NC_003076.8",
                                    "NW_003302551.1", "NW_003302554.1",
                                    "NW_003302553.1", "NW_003302552.1",
                                    "NW_003302551.1", "NW_003302550.1",
                                    "NW_003302549.1", "NW_003302548.1"))
g
## change title size
g + ggplot2::theme(title = element_text(size = 16))
## change axis title size (gene1.mid; gene2.mid)
g + ggplot2::theme(axis.title.x = element_text(size = 16),
                   axis.title.y = element_text(size = 16))
## change grid title size
g + theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16))
## change grid axis size and angle
g + theme(axis.text.x = element_text(size = 12, angle = 90))

## get help
#?plot_dagchainer

## -----------------------------------------------------------------------------
## load example Ka/Ks values (see above commands above how to obtain them)
data("ath_aly_ncbi_kaks", package="CRBHits")

## plot Ka/Ks results as histogram colored by Ka/Ks values
g <- plot_kaks(ath_aly_ncbi_kaks)

## plot Ka/Ks results as histogram filter for ka.min, ka.max, ks.min, ks.max
g.min_max <- plot_kaks(ath_aly_ncbi_kaks, ka.min = 0, ka.max = 1, ks.min = 0, ks.max = 1)

## select subset of chromosomes - needs gene position information
head(ARATHA.cds.longest.genepos)
head(ARALYR.cds.longest.genepos)
g.subset <- plot_kaks(ath_aly_ncbi_kaks,
                      gene.position.cds1 = ARATHA.cds.longest.genepos,
                      gene.position.cds2 = ARALYR.cds.longest.genepos,
                      select.chr = c("NC_003070.9", "NC_003071.7",
                        "NW_003302551.1", "NW_003302554.1"))

## plot Ka/Ks results and split by chromosome - needs gene position information

g.split <- plot_kaks(ath_aly_ncbi_kaks,
                     gene.position.cds1 = ARATHA.cds.longest.genepos,
                     gene.position.cds2 = ARALYR.cds.longest.genepos,
                     select.chr = c("NC_003070.9", "NC_003071.7",
                        "NW_003302551.1", "NW_003302554.1"),
                     splitByChr = TRUE)

## -----------------------------------------------------------------------------
## examples how to filter and mutate Ka/Ks results

## filter for Ks values < 2 on original Ka/Ks results
head(ath_aly_ncbi_kaks)
dim(ath_aly_ncbi_kaks)
dim(ath_aly_ncbi_kaks %>% dplyr::filter(ks < 2))

## filter for Ks values < 2 on plot object
head(g.split$g.kaks$data)
dim(g.split$g.kaks$data)
dim(g.split$g.kaks$data %>% dplyr::filter(ks < 2))

## filter for Ks values < 1 on plot object and plot
g.split$g.kaks$data %>% dplyr::filter(ks < 1) %>%
  ggplot2::ggplot() + ggplot2::geom_histogram(binwidth = 0.1, aes(x = ks))

## -----------------------------------------------------------------------------
## example comparing Homo sapiens and Pan troglodytes

## set URLs for Homo sapiens and Pan troglodytes from Ensembl

## set Ensembl URL
ensembl <- "ftp://ftp.ensembl.org/pub/release-101/fasta/"
## set Homo sapiens CDS URL
HOMSAP.cds.url <- paste0(ensembl, "homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz")
## set Pan troglodytes CDS URL
PANTRO.cds.url <- paste0(ensembl, "pan_troglodytes/cds/Pan_troglodytes.Pan_tro_3.0.cds.all.fa.gz")

## get Homo sapiens CDS
HOMSAP.cds.file <- tempfile()
#download.file(HOMSAP.cds.url, HOMSAP.cds.file, quiet = FALSE)
curl::curl_download(url = HOMSAP.cds.url,
                    handle = ftp_handle,
                    mode = "wb",
                    destfile = HOMSAP.cds.file)
HOMSAP.cds <- Biostrings::readDNAStringSet(HOMSAP.cds.file)

## get Pan troglodytes CDS
PANTRO.cds.file <- tempfile()
#download.file(PANTRO.cds.url, PANTRO.cds.file, quiet = FALSE)
curl::curl_download(url = PANTRO.cds.url,
                    handle = ftp_handle,
                    mode = "wb",
                    destfile = PANTRO.cds.file)
PANTRO.cds <- Biostrings::readDNAStringSet(PANTRO.cds.file)

## get longest isoform
HOMSAP.cds.longest <- isoform2longest(HOMSAP.cds, "ENSEMBL")
PANTRO.cds.longest <- isoform2longest(PANTRO.cds, "ENSEMBL")

## calculate CRBHit pairs for H. sapiens and P. troglodytes using 2 threads

## input from CDS obtained from Ensembl
## longest isoform selection
## query coverage >= 50%
## rost199 filter
HOMSAP_PANTRO_crbh <- cds2rbh(HOMSAP.cds,
                              PANTRO.cds,
                              qcov = 0.5, rost1999 = TRUE,
                              longest.isoform = TRUE,
                              isoform.source = "ENSEMBL",
                              threads = 2, plotCurve = TRUE,
                              lastpath = vignette.paths[1])

## get gene position for H. sapiens longest isoforms
HOMSAP.cds.longest.genepos <- cds2genepos(HOMSAP.cds.longest, source = "ENSEMBL")
## get gene position for P. troglodytes longest isoforms
PANTRO.cds.longest.genepos <- cds2genepos(PANTRO.cds.longest, source = "ENSEMBL")

## DAGchainer using gene base pair (gene bp start end)
## Note: change parameter to fit bp option
HOMSAP_PANTRO_crbh.dagchainer.bp <- rbh2dagchainer(rbhpairs = HOMSAP_PANTRO_crbh, 
                                      gene.position.cds1 = HOMSAP.cds.longest.genepos,
                                      gene.position.cds2 = PANTRO.cds.longest.genepos,
                                      type = "bp",
                                      gap_length = 10000,
                                      max_dist_allowed = 200000,
                                      dagchainerpath = vignette.paths[3])

## plot DAGchainer results for selected chromosome combinations
plot_dagchainer(HOMSAP_PANTRO_crbh.dagchainer.bp, select.chr = c("1","2","3","4","5","14"))

## -----------------------------------------------------------------------------
## load example Ka/Ks values (see above commands above how to obtain them)
data("hom_pan_ensembl_kaks", package="CRBHits")

## plot Ka/Ks results as histogram colored by Ka/Ks values
g <- plot_kaks(hom_pan_ensembl_kaks)

