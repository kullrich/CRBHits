#' @title cdsdir2orthofinder
#' @name cdsdir2orthofinder
#' @description This function calculates (conditional-)reciprocal best hit
#' (CRBHit) pairs for all possible comparison including self
#' comparison from a directory of CDS fasta files.
#' Sequence searches are performed with \bold{last}
#' \emph{Kiełbasa, SM et al. (2011)} [default]
#' or with \bold{mmseqs2}
#' \emph{Steinegger, M and Soeding, J (2017)}
#' or with \bold{diamond}
#' \emph{Buchfink, B et al. (2021)}.
#' @param dir directory containing CDS fasta files [mandatory]
#' @param file_ending define file ending to consider [default: *]
#' @param searchtool specify sequence search algorithm last, mmseqs2, diamond or
#' lambda3
#' [default: last]
#' @param lastpath specify the PATH to the last binaries
#' [default: /extdata/last-1595/bin/]
#' @param lastD last option D: query letters per random alignment
#' [default: 1e6]
#' @param lastm last option m: maximum initial matches per query position
#' [default: 10]
#' @param mmseqs2path specify the PATH to the mmseqs2 binaries
#' [default: NULL]
#' @param mmseqs2sensitivity specify the sensitivity option of mmseqs2
#' [default: 5.7]
#' @param mmseqs2maxseqs mmseqs2 option: Maximum results per query sequence
#' allowed to pass the prefilter
#' [default: 300]
#' @param diamondpath specify the PATH to the diamond binaries
#' [default: NULL]
#' @param diamondsensitivity specify the sensitivity option of diamond
#' [default: --sensitive]
#' @param diamondmaxtargetseqs specify the maximum number of target sequences
#' per query option of diamond
#' [default: 0]
#' @param lambda3path specify the PATH to the lambda3 binaries
#' [default: NULL]
#' @param lambda3sensitivity specify the sensitivity option of lambda3
#' [default: sensitive]
#' @param lambda3nummatches specify the number of matches per query option of
#' lambda3
#' [default: 25]
#' @param outpath specify the output PATH [default: /tmp]
#' @param dbpath specify the output PATH [default: /tmp]
#' @param crbh specify if conditional-reciprocal hit pairs should be retained
#' as secondary hits [default: TRUE]
#' @param keepSingleDirection specify if single direction secondary hit pairs
#' should be retained [default: FALSE]
#' @param eval evalue [default: 1e-3]
#' @param qcov query coverage [default: 0.0]
#' @param tcov target coverage [default: 0.0]
#' @param pident percent identity [default: 0.0]
#' @param alnlen alignment length [default: 0.0]
#' @param rost1999 specify if hit pairs should be filter by equation 2 of
#' Rost 1999 [default: FALSE]
#' @param filter specify additional custom filters as list to be applied on
#' hit pairs [default: NULL]
#' @param fit.type specify if mean or median should be used for fitting
#' [default: mean]
#' @param fit.varweight factor for fitting function to consider neighborhood
#' [default: 0.1]
#' @param fit.min specify minimum neighborhood alignment length
#' [default: 5]
#' @param threads number of parallel threads [default: 1]
#' @param remove specify if last result files should be removed [default: TRUE]
#' @param remove.db specify if last db files should be removed [default: TRUE]
#' @return List of three (crbh=FALSE)\cr
#' @importFrom Biostrings writeXStringSet
#' @importFrom graphics legend par points
#' @importFrom stats median splinefun
#' @importFrom utils read.table combn
#' @importFrom tidyr %>%
#' @seealso \code{\link[CRBHits]{aafile2rbh}}
#' @references Aubry S, Kelly S et al. (2014) Deep Evolutionary Comparison of
#' Gene Expression Identifies Parallel Recruitment of Trans-Factors in Two
#' Independent Origins of C4 Photosynthesis. \emph{PLOS Genetics},
#' \bold{10(6)} e1004365.
#' @references Kiełbasa, SM et al. (2011) Adaptive seeds tame genomic sequence
#' comparison. \emph{Genome research}, \bold{21(3)}, 487-493.
#' @references Rost B. (1999). Twilight zone of protein sequence alignments.
#' \emph{Protein Engineering}, \bold{12(2)}, 85-94.
#' @examples
#' ## compile last-1542 within CRBHits
#' CRBHits::make_last()
#' @export cdsdir2orthofinder
#' @author Kristian K Ullrich

cdsdir2orthofinder <- function(dir,
    file_ending="*",
    searchtool="last",
    lastpath=paste0(find.package("CRBHits"),
        "/extdata/last-1595/bin/"),
    lastD=1e6,
    lastm=10,
    mmseqs2path=NULL,
    mmseqs2sensitivity=5.7,
    mmseqs2maxseqs=300,
    diamondpath=NULL,
    diamondsensitivity="--sensitive",
    diamondmaxtargetseqs=0,
    lambda3path=NULL,
    lambda3sensitivity="sensitive",
    lambda3nummatches=25,
    outpath="/tmp",
    dbpath="/tmp",
    crbh=TRUE,
    keepSingleDirection=FALSE,
    eval=1e-3,
    qcov=0.0,
    tcov=0.0,
    pident=0.0,
    alnlen=0.0,
    rost1999=FALSE,
    filter=NULL,
    fit.type="mean",
    fit.varweight=0.1,
    fit.min=5,
    threads=1,
    remove=FALSE,
    remove.db=FALSE
    ){
    if(searchtool=="last"){
        if(!dir.exists(lastpath)){
            stop("Error: last PATH does not exist. Please specify correct
                PATH and/or look into package installation prerequisites.
                Try to use make_last() function.")
        }
        if(!file.exists(paste0(lastpath, "lastdb"))){
            stop("Error: lastdb binary does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites. Try to use make_last() function.")
        }
        if(!file.exists(paste0(lastpath, "lastal"))){
            stop("Error: lastal binary does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites. Try to use make_last() function.")
        }
    }
    if(searchtool=="mmseqs2"){
        if(!dir.exists(mmseqs2path)){
            stop("Error: mmseqs2 PATH does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
        if(!file.exists(paste0(mmseqs2path, "mmseqs"))){
            stop("Error: mmseqs2 binary does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
    }
    if(searchtool=="diamond"){
        if(!dir.exists(diamondpath)){
            stop("Error: diamond PATH does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
        if(!file.exists(paste0(diamondpath, "diamond"))){
            stop("Error: diamond binary does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
    }
    if(searchtool=="lambda3"){
        if(!dir.exists(lambda3path)){
            stop("Error: lambda3 PATH does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
        if(!file.exists(paste0(lambda3path, "lambda3"))){
            stop("Error: lambda3 binary does not exist. Please specify
                correct PATH and/or look into package installation
                prerequisites.")
        }
    }
    if(!dir.exists(outpath)){
        dir.create(outpath)
    }
    if(!dir.exists(dbpath)){
        dir.create(dbpath)
    }
    cds_files <- list.files(dir, file_ending)
    cds_species <- seq_along(cds_files)-1
    cds_species_files <- paste0("Species", cds_species, ".cds.fa")
    aa_species_files <- paste0("Species", cds_species, ".aa.fa")
    aa_db_files <- paste0("Species", cds_species, ".db")
    # write SpeciesIDs.txt
    sink(file.path(outpath, "SpeciesIDs.txt"))
    cat(
        apply(cbind(cds_species, cds_files), 1,
            function(x){paste0(x[1],": ", x[2])}), sep="\n"
    )
    sink(NULL)
    # extract sequence names, create cds_species_files
    # and perform selfblast
    for(i in seq_along(cds_files)){
        tmp_cds <- Biostrings::readDNAStringSet(
            file.path(dir, cds_files[i]))
        tmp_cds_names <- paste0(
            cds_species[i], "_", seq_along(tmp_cds)-1)
        cds_sequence_names <- data.frame(tmp_cds_names,
            stringr::word(names(tmp_cds)))
        # write SequenceIDs.txt
        write.table(cds_sequence_names, sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=TRUE,
            file=file.path(outpath, "SequenceIDs.txt"))
        names(tmp_cds) <- tmp_cds_names
        Biostrings::writeXStringSet(tmp_cds,
            file.path(outpath, cds_species_files[i]))
        tmp_aa <- MSA2dist::cds2aa(
            Biostrings::DNAStringSet(
                gsub("\\.", "-", tmp_cds)), shorten = TRUE)
        Biostrings::writeXStringSet(tmp_aa,
            file.path(outpath, aa_species_files[i]))
        tmp_aa_rbh <- aafile2rbh(
            aafile1 = file.path(outpath, aa_species_files[i]),
            aafile2 = file.path(outpath, aa_species_files[i]),
            dbfile1 = file.path(dbpath, aa_db_files[i]),
            dbfile2 = file.path(dbpath, aa_db_files[i]),
            searchtool=searchtool,
            lastpath=lastpath,
            lastD=lastD,
            lastm=lastm,
            mmseqs2path=mmseqs2path,
            mmseqs2sensitivity=mmseqs2sensitivity,
            mmseqs2maxseqs=mmseqs2maxseqs,
            diamondpath=diamondpath,
            diamondsensitivity=diamondsensitivity,
            diamondmaxtargetseqs=diamondmaxtargetseqs,
            lambda3path=lambda3path,
            lambda3sensitivity=lambda3sensitivity,
            lambda3nummatches=lambda3nummatches,
            outpath=outpath,
            crbh=crbh,
            keepSingleDirection=keepSingleDirection,
            eval=eval,
            qcov=qcov,
            tcov=tcov,
            pident=pident,
            alnlen=alnlen,
            rost1999=rost1999,
            filter=filter,
            plotCurve=FALSE,
            fit.type=fit.type,
            fit.varweight=fit.varweight,
            fit.min=fit.min,
            threads=threads,
            aafile2tmp=FALSE,
            remove=remove,
            remove.db=remove.db)
        write.table(tmp_aa_rbh$crbh1[,1:12], sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE,
            file=file.path(outpath, paste0("Blast", cds_species[i], "_",
            cds_species[i], ".txt")))
    }
    # combn
    to_calc <- t(combn(cds_species, 2))
    apply(to_calc, 1, function(x){
        tmp_aa_rbh <- aafile2rbh(
            aafile1 = file.path(outpath, aa_species_files[x[1]+1]),
            aafile2 = file.path(outpath, aa_species_files[x[2]+1]),
            dbfile1 = file.path(dbpath, aa_db_files[x[1]+1]),
            dbfile2 = file.path(dbpath, aa_db_files[x[2]+1]),
            searchtool=searchtool,
            lastpath=lastpath,
            lastD=lastD,
            lastm=lastm,
            mmseqs2path=mmseqs2path,
            mmseqs2sensitivity=mmseqs2sensitivity,
            mmseqs2maxseqs=mmseqs2maxseqs,
            diamondpath=diamondpath,
            diamondsensitivity=diamondsensitivity,
            diamondmaxtargetseqs=diamondmaxtargetseqs,
            lambda3path=lambda3path,
            lambda3sensitivity=lambda3sensitivity,
            lambda3nummatches=lambda3nummatches,
            outpath=outpath,
            crbh=crbh,
            keepSingleDirection=keepSingleDirection,
            eval=eval,
            qcov=qcov,
            tcov=tcov,
            pident=pident,
            alnlen=alnlen,
            rost1999=rost1999,
            filter=filter,
            plotCurve=FALSE,
            fit.type=fit.type,
            fit.varweight=fit.varweight,
            fit.min=fit.min,
            threads=threads,
            aafile2tmp=FALSE,
            remove=remove,
            remove.db=remove.db)
        write.table(tmp_aa_rbh$crbh1[,1:12], sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE,
            file=file.path(outpath, paste0("Blast", x[1], "_", x[2], ".txt")))
        write.table(tmp_aa_rbh$crbh2[,1:12], sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE,
            file=file.path(outpath, paste0("Blast", x[2], "_", x[1], ".txt")))
    })
}
