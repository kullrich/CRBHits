[![Build Status](https://travis-ci.com/kullrich/CRBHits.svg?token=CdAqxcdzRKt5DsRtbcmR&branch=master)](https://travis-ci.com/kullrich/CRBHits)
[![pipeline status](https://gitlab.gwdg.de/mpievolbio-it/crbhits/badges/master/pipeline.svg)](https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/commits/master)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![MITlicense](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

CRBHits - From Conditional Reciprocal Best Hits to Codon Alignments and Ka/Ks in R
=========

R package source code: [https://gitlab.gwdg.de/mpievolbio-it/crbhits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

R package pages: [https://mpievolbio-it.pages.gwdg.de/crbhits/](https://mpievolbio-it.pages.gwdg.de/crbhits/)

R package issues: [https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues](https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues)

# CRBHits - Description

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is a reimplementation of the Conditional Reciprocal Best Hit (CRBH) algorithm [crb-blast](https://github.com/cboursnell/crb-blast) in R. It covers all necessary steps from CRBHit pair calculation to Codon Alignments and Ka/Ks.

The CRBH algorithm was introduced by [Aubry S, Kelly S et al. (2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004365) and ported to python [shmlast](https://pypi.org/project/shmlast/) ([Scott C. 2017](https://joss.theoj.org/papers/10.21105/joss.00142)) which benefits from the blast-like sequence search software [LAST](http://last.cbrc.jp/) ([Kiełbasa SM et al. 2011](https://genome.cshlp.org/content/21/3/487.long)).

Like [shmlast](https://pypi.org/project/shmlast/), [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) plots the fitted model of the CRBHit evalue based algorithm. In addition users can filter the CRBHit pairs prior fitting for other criteria like evalue, protein identity and/or the twilight zone of protein sequence alignments according to [Rost B. (1999)](https://academic.oup.com/peds/article/12/2/85/1550637).

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) only takes coding nucleotide sequences as the query and target inputs, since the secondary aim of CRBHits is to calculate synonymous and non-synonymous substitutions with the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) ([Sharif D, Lobry JR. 2007](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)) or the external tool [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) ([Wang D, Zhang Y et al. 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5054116/)).

This is in contrast to [crb-blast](https://github.com/cboursnell/crb-blast), which can take proteins or nucleotides as the target sequences and in contrast to [shmlast](https://pypi.org/project/shmlast/), which uses nucleotides as queries and proteins as target sequences.

The resulting CRBHit pairs can be used to obtain pairwise codon alignments, which are further used to calculate synonymous and nonsynonymous substitutions using parallelization.

The Ka/Ks (also sometimes denoted as dN/dS) values can be obtained either via the codon model of [Li WH. (1999)](https://www.ncbi.nlm.nih.gov/pubmed/8433381) as implemented in the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) or the model of [Yang Z and Nielson R. (2000)](https://www.ncbi.nlm.nih.gov/pubmed/10666704) as implemented in [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download).

The following two images show the two main functions of the package `cds2rbh()` and `rbh2kaks()`, which are described in more detail in the package vignettes.

<img src="./vignettes/cds2rbhoverview.png" alt="Figure: Overview of the `cds2rbh()` function" width="200">

<img src="./vignettes/rbh2kaksoverview.png" alt="Figure: Overview of the `rbh2kaks()` function" width="200">

# Installation

see also here for the R package pages [https://mpievolbio-it.pages.gwdg.de/crbhits/](https://mpievolbio-it.pages.gwdg.de/crbhits/)

## R specific installation prerequisites

### install packages from [cran](https://cran.r-project.org/web/packages/index.html)

In most cases you need to first install the following system-wide packages to be able to compile the R dependencies.

Ubuntu/Debian

```
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev libglu1-mesa-dev libgit2-dev
#pkgdown dependencies - pkgdown is used to build R package pages
#sudo apt-get install libssh2-1-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

CentOS

```
sudo yum install libcurl-devel openssl-devel libxml2-devel mesa-libGLU-devel libgit2-devel
#pkgdown dependencies - pkgdown is used to build R package pages
#sudo yum install libssh2-devel fontconfig-devel harfbuzz-devel fribidi-devel
```

- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
- [testthat](https://cran.r-project.org/web/packages/testthat/index.html)
- [curl](https://cran.r-project.org/web/packages/curl/index.html)
- [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](https://cran.r-project.org/web/packages/doMC/index.html)
- [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
- [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
- [mclust](https://cran.r-project.org/web/packages/mclust/index.html)
- [feature](https://cran.r-project.org/web/packages/feature/index.html)

```
install.packages("devtools")
install.packages("testthat")
install.packages("curl")
install.packages("seqinr")
install.packages("foreach")
install.packages("doMC")
install.packages("tidyverse")
install.packages("gridExtra")
install.packages("mclust")
install.packages("feature")
```

### install packages from [Bioconductor](https://www.bioconductor.org/)

- [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### install [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

```
library(devtools)
install_gitlab("mpievolbio-it/crbhits", host = "https://gitlab.gwdg.de",
build_vignettes = FALSE, dependencies = FALSE)
#install_github("kullrich/CRBHits", build_vignettes = FALSE, dependencies = FALSE)
```

## External tools installation prerequisites

The source code for the prerequisites (LAST, KaKs_Calculator2.0, DAGchainer) are forked within [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). 

- [LAST](http://last.cbrc.jp/) [http://last.cbrc.jp/last-1133.zip](http://last.cbrc.jp/last-1133.zip)

To compile the forked version of [LAST](http://last.cbrc.jp/) within the `CRBHits` R package directory try to use the function `make.last()`:

```
## see more installation information here, if make fails
## last-install-help: http://last.cbrc.jp/doc/last.html
make.last()
```

To compile [LAST](http://last.cbrc.jp/) yourself on Linux/Unix/macOS into another folder:

```
## create and change into the directory to install LAST
## e.g.
mkdir /tmp/last
cd /tmp/last
## donwload last-1133
curl -O http://last.cbrc.jp/last-1133.zip
unzip last-1133.zip
cd last-1133
## compile LAST
make
```

- [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)

To compile the forked version of [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) within the `CRBHits` R package directory try to use the function `make.KaKs_Calculator2()`:

```
## compile KaKs_Calculator2
make.KaKs_Calculator2()
```

__Note:__ Due to some changes in the latest **g++** compilers the source code was altered to meet this changes, which are directly incorporated into the `KaKs_Calculator2.0.tar.gz` that is distributed with [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). It is recommended to compile from this file (see below):

```
## create and change into the directory to install KaKs_Calculator2
## e.g.
mkdir /tmp/KaKs_Calculator2
cd /tmp/KaKs_Calculator2
## donwload KaKs_Calculator2
curl -O https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/raw/master/inst/extdata/KaKs_Calculator2.0.tar.gz
tar -xvf KaKs_Calculator2.0.tar.gz
cd KaKs_Calculator2.0/src
## compile KaKs_Calculator2
make clean
make
```

- [DAGchainer](http://dagchainer.sourceforge.net/)

To compile the forked version of [DAGchainer](http://dagchainer.sourceforge.net/) within the `CRBHits` R package directory try to use the function `make.dagchainer()`:

```
## compile DAGchainer
make.dagchainer()
```

__Note:__ Due to some changes in the latest **g++** compilers the source code was altered to meet this changes, which are directly incorporated into the `dagchainer.zip` that is distributed with [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). It is recommended to compile from this file (see below):

```
## create and change into the directory to install DAGchainer
## e.g.
mkdir /tmp/DAGchainer
cd /tmp/DAGchainer
## donwload DAGchainer
curl -O https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/raw/master/inst/extdata/dagchainer.zip
unzip dagchainer.zip
cd dagchainer
## compile DAGchainer
make
```

If you would like to install the latest version of the tools, you need to download the source code and compile again.

If you would like to use your own compiled versions of `LAST`, `KaKs_Calculator2.0` and `DAGchainer` you need to set the correct `@param` in the corresponding functions of [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits).

```
## example how to use own compiled versions of LAST, KaKs_Calculator2.0 and DAGchainer
my.lastpath <- "/tmp/last/last-1133/src"
my.kakspath <- "/tmp/KaKs_Calculator2/KaKs_Calculator2.0/src"
my.dagchainerpath <- "/tmp/dagcahiner"

?cds2rbh
cds2rbh(., ., lastpath = my.lastpath)

?rbh2kaks
rbh2kaks(., ., model = "YN", kakscalcpath = my.kakspaths)

?rbh2dagchainer
rbh2dagchainer(., ., dagchainerpath = my.dagchainerpath)
```

## Vignettes

These vignettes introduce  [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

- [CRBHits Basic Vignette](https://mpievolbio-it.pages.gwdg.de/crbhits/articles/V01CRBHitsBasicVignette.html) - Basic Usage of CRBHits
    - includes CRBHit pair calculation
    - includes CRBHit pair filtering
    - includes Longest Isoform selection
    - includes Codon alignments
    - includes Ka/Ks calculations
- [KaKs Vignette](https://mpievolbio-it.pages.gwdg.de/crbhits/articles/V02KaKsVignette.html) - KaKs Calculations between two species and subsequent data filter steps
    - includes *Arabidopsis thalina* and *Arabidopsis lyrata* CRBHit pair calculation
    - includes *Homo sapiens* and *Pan troglodytes* CRBHit pair calculation
    - inlcudes Longest Isoform selection
    - includes Gene/Isoform chromosomal position extraction
    - includes Tandem Duplicate Assignment
    - includes Synteny Assignment
    - includes Ka/Ks colored Dot-Plot

## Quick-guide

```
library(CRBHits)
## compile last-1133
make.last()
## conditional reciprocal best hits (CRBHit pairs)
data("ath", package="CRBHits")
data("aly", package="CRBHits")
ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE)
summary(ath_aly_crbh)

?cds2rbh

# #kaks calculation - subset model "Li"
ath_aly_crbh$crbh.pairs <- head(ath_aly_crbh$crbh.pairs, 20)
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh,
                              ath, aly, model = "Li")
head(ath_aly_crbh.kaks)
?rbh2kaks

## plot kaks
g.kaks <- plot.kaks(ath_aly_crbh.kaks)
?plot.kaks

## kaks calculation - subset model "YN"
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh,
                              ath, aly, model = "YN")

## kaks calculation using multiple threads
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh, ath, aly, threads = 2)
head(ath_aly_crbh.kaks)

## example how to use own compiled versions of LAST
my.lastpath <- "/tmp/last/last-1133/src"
ath_aly_crbh <- cds2rbh(ath, aly, plotCurve = TRUE,
                        lastpath = my.lastpath)
?cds2rbh

## example how to use own compiled versions of KaKs_Calculator2.0
my.kakspath <- "/tmp/KaKs_Calculator2/KaKs_Calculator2.0/src"
ath_aly_crbh.kaks <- rbh2kaks(ath_aly_crbh,
                              ath, aly, model = "YN",
                              kakscalcpath = my.kakspath)
?rbh2kaks
```

## License

MIT (see LICENSE)

The [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) package includes source code that has been published under following licenses:

### last-1133.zip

GNU General Public License Version 3, 29 June 2007 [GPLv3](https://www.gnu.org/licenses/gpl-3.0.de.html)

### KaKs_Calculator2.0.tar.gz

The toolkit is freely available (licensed under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.de.html)) online at [https://sourceforge.net/projects/kakscalculator2/](https://sourceforge.net/projects/kakscalculator2/).

### DAGchainer

All of this software is open source, free of license restrictions. 

## Contributing Code

If you would like to contribute to CRBHits, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the CRBHits developer team agrees that it’s a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your contribution with package development as a whole, [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [main repo](https://gitlab.gwdg.de/mpievolbio-it/crbhits.git), branch off a [feature branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) from `master`, [commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project) and [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) your changes to your fork and submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) for `CRBHits:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please report any errors or requests regarding [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) to Kristian Ullrich (ullrich@evolbio.mpg.de)

or use the issue tracker at [https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues](https://gitlab.gwdg.de/mpievolbio-it/crbhits/issues)

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

## References

Aubry S., Kelly S., Kümpers B. M., Smith-Unna R. D., and Hibberd J. M. (2014). **Deep evolutionary comparison of gene expression identifies parallel recruitment of trans-factors in two independent origins of C4 photosynthesis.** *PLoS Genetics*, **10(6)**. [https://doi.org/10.1371/journal.pgen.1004365](https://doi.org/10.1371/journal.pgen.1004365)

Charif D., and Lobry J. R. (2007). **SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.** *In Structural approaches to sequence evolution* (pp. 207-232). Springer, Berlin, Heidelberg. [https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)

Duong T., and Wand M. (2015). **feature: Local Inferential Feature Significance for Multivariate Kernel Density Estimation.** *R package version 1.2.13*. [https://cran.r-project.org/web/packages/feature/](https://cran.r-project.org/web/packages/feature/)

Haas B. J., Delcher A. L., Wortman J. R., and Salzberg S. L. (2004). **DAGchainer: a tool for mining segmental genome duplications and synteny.** *Bioinformatics*, **20(18)**, 3643-3646. [https://doi.org/10.1093/bioinformatics/bth397](https://doi.org/10.1093/bioinformatics/bth397)

Haug-Baltzell A., Stephens S. A., Davey S., Scheidegger C. E., Lyons E. (2017). **SynMap2 and SynMap3D: web-based wholge-genome synteny browsers.** *Bioinformatics*, **33(14)**, 2197-2198. [https://academic.oup.com/bioinformatics/article/33/14/2197/3072872](https://academic.oup.com/bioinformatics/article/33/14/2197/3072872)

Kiełbasa S. M., Wan R., Sato K., Horton P., and Frith M. C. (2011). **Adaptive seeds tame genomic sequence comparison.** *Genome Research*, **21(3)**, 487-493. [https://doi.org/10.1101/gr.113985.110](https://doi.org/10.1101/gr.113985.110)

Kimura M. (1977). **Preponderance of synonymous changes as evidence for the neutral theory of molecular evolution.** *Nature*, **267**, 275-276.

Li W. H. (1993). **Unbiased estimation of the rates of synonymous and nonsynonymous substitution.** *Journal of Molecular Evolution*, **36(1)**, 96-99. [https://doi.org/10.1007/bf02407308](https://doi.org/10.1007/bf02407308)

Microsoft, and Weston S. (2020). **foreach: Provides Foreach Looping Construct.** *R package version*, 1.5.1. [foreach](https://cran.r-project.org/web/packages/foreach/index.html)

Mugal C. F., Wolf J. B. W., Kaj I. (2014). **Why Time Matters: Codon Evolution and the Temproal Dynamics of dN/dS.** *Molecular Biology and Evolution*, **31(1)**, 212-231.

Ooms J. (2019). **curl: A Modern and Flexible Web Client for R.** *R package version*, 4.3. [curl](https://cran.r-project.org/web/packages/curl/index.html)

Pagès H., Aboyoun P., Gentleman R., and DebRoy S. (2017). **Biostrings: Efficient manipulation of biological strings.** *R package version*, 2.56.0. [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html)

Revolution Analytics, and Weston S. (2020). **doMC: Foreach Parallel Adaptor for 'parallel'.** *R package version*, 1.3.7. [doMC](https://cran.r-project.org/web/packages/doMC/index.html)

Rost B. (1999). **Twilight zone of protein sequence alignments.** *Protein Engineering*, **12(2)**, 85-94. [https://doi.org/10.1093/protein/12.2.85](https://doi.org/10.1093/protein/12.2.85)

Scott C. (2017). **shmlast: an improved implementation of conditional reciprocal best hits with LAST and Python.** *Journal of Open Source Software*, **2(9)**, 142. [https://joss.theoj.org/papers/10.21105/joss.00142](https://joss.theoj.org/papers/10.21105/joss.00142)

Scrucca L., Fop M., Murphy T. B., and Raftery A. E. (2016) **mclust 5: clustering, classification and density estimation using Gaussian finite mixture models.** *The R Journal*, 8(1), 289-317. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/)

Wickham H. (2011). **testthat: Get Started with Testing.** *The R Journal*, 3(1), 5. [testthat](https://cran.r-project.org/web/packages/testthat/index.html)

Wickham H. (2019). **stringr: Simple, Consistent Wrappers for Common String Operations.** *R package version*, 1.4.0. [stringr](https://cran.r-project.org/web/packages/stringr/index.html)

Wickham H. (2020). **tidyr: Tidy Messy Data.** *R package version*, 1.1.2. [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)

Wickham H., Hester J., and Chang W. (2020). **devtools: Tools to make Developing R Packages Easier.** *R package version*, (2.3.2). [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

Wickham H., François R., Henry L., and Müller K. (2020). **dplyr: A Grammar of Data Manipulation.** *R package version*, 1.0.2. [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)

Yang Z., and Nielsen R. (2000). **Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models.** *Molecular Biology and Evolution*, **17(1)**, 32-43. [https://doi.org/10.1093/oxfordjournals.molbev.a026236](https://doi.org/10.1093/oxfordjournals.molbev.a026236)
