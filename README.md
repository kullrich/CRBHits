[![Build Status](https://travis-ci.com/kullrich/CRBHits.svg?token=CdAqxcdzRKt5DsRtbcmR&branch=master)](https://travis-ci.com/kullrich/CRBHits)
[![pipeline status](https://gitlab.gwdg.de/mpievolbio-it/crbhits/badges/master/pipeline.svg)](https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/commits/master)

CRBHits - Conditional Reciprocal Best Hits in R
=========

R package source code: [https://gitlab.gwdg.de/mpievolbio-it/crbhits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

R package pages: [https://mpievolbio-it.pages.gwdg.de/crbhits/](https://mpievolbio-it.pages.gwdg.de/crbhits/)

R packges issues: [https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues](https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues)

# CRBHits - Description

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is a reimplementation of the Conditional Reciprocal Best Hit algorithm [crb-blast](https://github.com/cboursnell/crb-blast) in R.

The algorithm was introduced by [Aubry S, Kelly S et al. (2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004365) and ported to python [shmlast](https://pypi.org/project/shmlast/) ([Scott C. 2017](https://joss.theoj.org/papers/10.21105/joss.00142)) which benefits from the blast-like sequence search software [LAST](http://last.cbrc.jp/) ([Kiełbasa SM et al. 2011](https://genome.cshlp.org/content/21/3/487.long)).

Like [shmlast](https://pypi.org/project/shmlast/), [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) plots the fitted model of the conditional reciprocal best hit evalue based algorithm. In addition users can filter the hit pairs prior fitting for other cirteria like evalue, protein identity and/or the twilight zone of protein sequence alignments according to [Rost B. (1999)](https://academic.oup.com/peds/article/12/2/85/1550637).

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) only takes coding nucleotide sequences as the query and target inputs, since the secondary aim of CRBHits is to calculate synonymous and non-synonymous substitutions with the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) ([Sharif D, Lobry JR. 2007](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)) or the external tool [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) ([Wang D, Zhang Y et al. 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5054116/)). This is in contrast to [crb-blast](https://github.com/cboursnell/crb-blast), which can take proteins or nucleotides as the target sequences and in contrast to [shmlast](https://pypi.org/project/shmlast/), which uses nucleotides as queries and proteins as target sequences.

The resultings conditional reciprocal best hit pair matrix can be used to obtain pairwise codon alignments, which are further used to calculate synonymous and nonsynonymous substitutions using parallelization. The ka/ks (dN/dS) values can be obtained either via the codon model of [Li WH. (1999)](https://www.ncbi.nlm.nih.gov/pubmed/8433381) implemented in the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) or the model [Yang Z and Nielson R. (2000)](https://www.ncbi.nlm.nih.gov/pubmed/10666704) implemented in [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download).

# Installation

see also here for the R package pages [https://mpievolbio-it.pages.gwdg.de/crbhits/](https://mpievolbio-it.pages.gwdg.de/crbhits/)

## R specific installation prerequisites

### install packages from [cran](https://cran.r-project.org/web/packages/index.html)

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
- [testthat](https://cran.r-project.org/web/packages/testthat/index.html)
- [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](https://cran.r-project.org/web/packages/doMC/index.html)
- [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html)
- [mclust](https://cran.r-project.org/web/packages/mclust/index.html)
- [feature](https://cran.r-project.org/web/packages/feature/index.html)

```
install.packages("devtools")
install.packages("testthat")
install.packages("seqinr")
install.packages("foreach")
install.packages("doMC")
install.packages("magrittr")
install.packages("mclust")
install.packages("feature")
```

### install packages from [Bioconductor](https://www.bioconductor.org/)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### install [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

```
library(devtools)
install_gitlab("mpievolbio-it/crbhits", host = "https://gitlab.gwdg.de", build_vignettes = TRUE, dependencies = FALSE)
#install_github("kullrich/CRBHits", build_vignettes = TRUE, dependencies = FALSE)
```

## External tools installation prerequisites

The source code for both prerequisites (LAST, KaKs_Calculator2.0) are forked within [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). 

- [LAST](http://last.cbrc.jp/) [http://last.cbrc.jp/last-1060.zip](http://last.cbrc.jp/last-1060.zip)

To compile the forked version of [LAST](http://last.cbrc.jp/) within the `CRBHits` R package directory try to use the function `make.last()`:

```
## see more installation information here, if make fails
## last-install-help: http://last.cbrc.jp/doc/last.html
make.last()
```

To compile [LAST](http://last.cbrc.jp/) yourself on Linux/Unix/macOS into another folder:

```
#create and change into the directory to install LAST
#e.g.
mkdir /tmp/last
cd /tmp/last
#donwload last-1060
curl -O http://last.cbrc.jp/last-1060.zip
unzip last-1060.zip
cd last-1060
#compile LAST
make
```

To compile the forked version of [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) within the `CRBHits` R package directory try to use the function `make.KaKs_Calculator2()`:

__Note:__ The original files can be downloaded here:

- [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)

```
make.KaKs_Calculator2()
```

__Note:__ Due to some changes in the latest **g++** compilers the source code was altered to meet this changes, which are directly incorporated into the `KaKs_Calculator2.0.tar.gz` that is distributed with [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). It is recommended to compile from this file (see below):

```
#create and change into the directory to install KaKs_Calculator2
#e.g.
mkdir /tmp/KaKs_Calculator2
cd /tmp/KaKs_Calculator2
#donwload KaKs_Calculator2
curl -O https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/raw/master/inst/extdata/KaKs_Calculator2.0.tar.gz
tar -xvf KaKs_Calculator2.0.tar.gz
cd KaKs_Calculator2.0/src
#compile KaKs_Calculator2
make clean
make
```

If you would like to install the latest version of both tools, you need to download the source code and compile again.

If you would like to use your own compiled versions of `LAST` and `KaKs_Calculator2.0` you need to set the correct `@param` in the corresponding functions of [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits).

```
#example how to use own compiled versions of LAST and KaKs_Calculator2.0
my.lastpath <- "/tmp/last/last-1060/src"
my.kakspath <- "/tmp/KaKs_Calculator2/KaKs_Calculator2.0/src"
?cds2rbh
cds2rbh(., ., lastpath = my.lastpath)
?rbh2kaks
rbh2kaks(., ., model = "YN", kakscalcpath = my.kakspaths)
```

## Vignettes

These vignettes introduce  [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

- [CRBHits Basic Vignette](https://mpievolbio-it.pages.gwdg.de/crbhits/articles/CRBHitsBasicVignette.html) - Basic Usage of CRBHits
- [dNdS Clustering Vignette](https://mpievolbio-it.pages.gwdg.de/crbhits/articles/dNdSClusteringVignette.html) - dNdS Calculations between two species (*A. thalina* and *A. lyrata*)
- [Paranome-based Whole Genome Duplication Vignette](https://mpievolbio-it.pages.gwdg.de/crbhits/articles/WGDVignette.html) - Selfblast and subsequent WGD prediction (*A. thaliana*)

## Quick-guide

```
library(CRBHits)
#compile last-1060
make.last()
#conditional-reciprocal best hits
data("ath", package="CRBHits")
data("aly", package="CRBHits")
ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE)
summary(ath_aly_crbh)
?cds2rbh

#kaks calculation
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh$crbh.pairs[1:20, ], ath, aly)
head(ath_aly_crbh.kaks)
?rbh2kaks

#kaks calculation using multiple threads
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh$crbh.pairs[1:20, ], ath, aly, threads = 4)
head(ath_aly_crbh.kaks)

#example how to use own compiled versions of LAST
my.lastpath <- "/tmp/last/last-1060/src"
?cds2rbh
ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE, lastpath = my.lastpath)

#example how to use own compiled versions of KaKs_Calculator2.0
my.kakspath <- "/tmp/KaKs_Calculator2/KaKs_Calculator2.0/src"
?rbh2kaks
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh$crbh.pairs[1:20, ],
                              ath, aly, model = "YN", kakscalcpath = my.kakspath)
```

## Todo
- kaks-clustering

## License

MIT (see LICENSE)

The [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) package includes source code that has been published under following licenses:

### last-1060.zip

GNU General Public License Version 3, 29 June 2007 [GPLv3](https://www.gnu.org/licenses/gpl-3.0.de.html)

### KaKs_Calculator2.0.tar.gz

The toolkit is freely available (licensed under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.de.html)) online at [https://sourceforge.net/projects/kakscalculator2/](https://sourceforge.net/projects/kakscalculator2/).

## Bug reports

Please report any errors or requests regarding [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) to Kristian Ullrich (ullrich@evolbio.mpg.de)

or use the issue tracker at [https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues](https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues)

## References

Aubry S., Kelly S., Kümpers B. M., Smith-Unna R. D., and Hibberd J. M. (2014). **Deep evolutionary comparison of gene expression identifies parallel recruitment of trans-factors in two independent origins of C4 photosynthesis.** *PLoS genetics*, **10(6)**. [https://doi.org/10.1371/journal.pgen.1004365](https://doi.org/10.1371/journal.pgen.1004365)

Charif D., and Lobry J. R. (2007). **SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.** *In Structural approaches to sequence evolution* (pp. 207-232). Springer, Berlin, Heidelberg. [https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)

Kiełbasa S. M., Wan R., Sato K., Horton P., and Frith M. C. (2011). **Adaptive seeds tame genomic sequence comparison.** *Genome research*, 21(3), 487-493. [https://doi.org/10.1101/gr.113985.110](https://doi.org/10.1101/gr.113985.110)

Li W. H. (1993). **Unbiased estimation of the rates of synonymous and nonsynonymous substitution.** *Journal of molecular evolution*, 36(1), 96-99. [https://doi.org/10.1007/bf02407308](https://doi.org/10.1007/bf02407308)

Pagès H., Aboyoun P., Gentleman R., and DebRoy S. (2017). **Biostrings: Efficient manipulation of biological strings.** *R package version*, 2(0).

Rost B. (1999). **Twilight zone of protein sequence alignments.** *Protein engineering*, 12(2), 85-94. [https://doi.org/10.1093/protein/12.2.85](https://doi.org/10.1093/protein/12.2.85)

Scott C. (2017). **shmlast: an improved implementation of conditional reciprocal best hits with LAST and Python.** *Journal of Open Source Software*, 2(9), 142. [https://joss.theoj.org/papers/10.21105/joss.00142](https://joss.theoj.org/papers/10.21105/joss.00142)

Scrucca L., Fop M., Murphy T. B., and Raftery A. E. (2016) **mclust 5: clustering, classification and density estimation using Gaussian finite mixture models.** *The R Journal*, 8(1), 289-317. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/)

Duong T., and Wand M. (2015). **feature: Local Inferential Feature Significance for Multivariate Kernel Density Estimation.** *R package version 1.2.13*. [https://cran.r-project.org/web/packages/feature/](https://cran.r-project.org/web/packages/feature/) 

Yang Z., and Nielsen R. (2000). **Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models.** *Molecular biology and evolution*, 17(1), 32-43. [https://doi.org/10.1093/oxfordjournals.molbev.a026236](https://doi.org/10.1093/oxfordjournals.molbev.a026236)
