CRBHits - Conditional Reciprocal Best Hit in R
=========
# CRBHits - Description

[CRBHits](https://github.com/kullrich/CRBHits) is a reimplementation of the Conditional Reciprocal Best Hit algorithm [crb-blast](https://github.com/cboursnell/crb-blast) in R.

The algorithm was introduced by [Aubry S, Kelly S et al 2014](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004365) and ported to python ([shmlast](https://pypi.org/project/shmlast/)) by [Scott C 2017](https://joss.theoj.org/papers/10.21105/joss.00142) which benefits from an alternative blast software, namely [LAST](http://last.cbrc.jp/) ([Kiełbasa SM et al 2011](https://genome.cshlp.org/content/21/3/487.long)).

Like [shmlast](https://pypi.org/project/shmlast/), [CRBHits](https://github.com/kullrich/CRBHits) plots the fitted model of the conditional reciprocal best hit evalue based algorithm. In addition users can filter the hit pairs prior fitting for other cirteria like evalue, protein identity and/or the twilight zone of protein sequence alignments according to [Rost 1999](https://academic.oup.com/peds/article/12/2/85/1550637).

[CRBHits](https://github.com/kullrich/CRBHits) only takes coding nucleotide sequences as the query and target inputs, since the primary aim of CRBHits is to calculate synonymous and non-synonymous substitutions with the R-package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) ([Sharif D, Lobry JR 2007](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)) or the external tool [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) ([Wang D, Zhang Y et al 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5054116/)). This is in contrast to [crb-blast](https://github.com/cboursnell/crb-blast), which can take proteins or nucleotides as the target sequences and in contrast to [shmlast](https://pypi.org/project/shmlast/), which 

# Installation

## R specific installation prerequisites

### install packages from cran [CRAN](https://cran.r-project.org/web/packages/index.html).

In most cases you need to first install the following system-wide packages to be able to compile the R dependencies.

Ubuntu/Debian

```
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

CentOS

```
sudo yum install libcurl-devel openssl-devel libxml2-devel
```

- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
- [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](https://cran.r-project.org/web/packages/doMC/index.html)
- [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html)

```
install.packages("devtools")
install.packages("seqinr")
install.packages("foreach")
install.packages("doMC")
install.packages("magrittr")
```

### install packages from Bioconductor [Bioconductor](https://www.bioconductor.org/).

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### install [CRBHits](https://github.com/kullrich/CRBHits)

```
library(devtools)
install_github("kullrich/CRBHits", build_vignettes = TRUE, dependencies = FALSE)
#install_git("https://gitlab.gwdg.de/mpievolbio-it/crbhits/crbhits.git", build_vignettes = TRUE, dependencies = FALSE)
```

## External tools installation prerequisites

The source code for both prerequisites (LAST, KaKs_Calculator2.0) are forked within [CRBHits](). 

- [LAST](http://last.cbrc.jp/)
```
#system("wget -O http://last.cbrc.jp/last-1060.zip -P ", paste0(find.package(CRBH), "/extdata/"))
setwd(paste0(find.package("CRBHits"), "/extdata/"))
system(paste0("unzip last-1060.zip"))
setwd(paste0(find.package("CRBHits"), "/extdata/last-1060/"))
# see more installation iformation here, if make fails
# last-install-help: http://last.cbrc.jp/doc/last.html
system("make")
```

- [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)

```
setwd(paste0(find.package("CRBHits"), "/extdata/"))
system(paste0("tar -xvf KaKs_Calculator2.0.tar.gz"))
setwd(paste0(find.package("CRBHits"), "/extdata/KaKs_Calculator2.0/src/"))
system("make clean")
system("make")
```

If you would like to install the latest version of both tools, you need to download the source code and compile again. 

## Vignettes

These vignettes introduce `CRBHits`

- [CRBHits Vignette](https://github.com/kullrich/CRBHits/tree/master/vignettes/CRBHitsVignette.Rmd)

## Quick-guide

```
library(CRBHits)
#conditional-reciprocal best hits
data("ath", package="CRBHits")
data("aly", package="CRBHits")
ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE)

#kaks calculation
ath_aly_crbh.kaks <- rbh2kaks(ath, aly, ath_aly_crbh)

#kaks calculation using multiple threads
ath_aly_crbh.kaks <- rbh2kaks(ath, aly, ath_aly_crbh, threads = 4)
```

## Todo
- 

## License

MIT (see LICENSE)

## Bug reports

Please report any errors or requests regarding `CRBHits` to Kristian Ullrich (ullrich@evolbio.mpg.de)
