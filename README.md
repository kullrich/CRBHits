CRBHits - Conditional Reciprocal Best Hits in R
=========
# CRBHits - Description

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is a reimplementation of the Conditional Reciprocal Best Hit algorithm [crb-blast](https://github.com/cboursnell/crb-blast) in R.

The algorithm was introduced by [Aubry S, Kelly S et al. (2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004365) and ported to python [shmlast](https://pypi.org/project/shmlast/) ([Scott C. 2017](https://joss.theoj.org/papers/10.21105/joss.00142)) which benefits from the blast-like sequence search software [LAST](http://last.cbrc.jp/) ([Kiełbasa SM et al. 2011](https://genome.cshlp.org/content/21/3/487.long)).

Like [shmlast](https://pypi.org/project/shmlast/), [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) plots the fitted model of the conditional reciprocal best hit evalue based algorithm. In addition users can filter the hit pairs prior fitting for other cirteria like evalue, protein identity and/or the twilight zone of protein sequence alignments according to [Rost B. (1999)](https://academic.oup.com/peds/article/12/2/85/1550637).

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) only takes coding nucleotide sequences as the query and target inputs, since the primary aim of CRBHits is to calculate synonymous and non-synonymous substitutions with the R-package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) ([Sharif D, Lobry JR. 2007](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)) or the external tool [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) ([Wang D, Zhang Y et al. 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5054116/)). This is in contrast to [crb-blast](https://github.com/cboursnell/crb-blast), which can take proteins or nucleotides as the target sequences and in contrast to [shmlast](https://pypi.org/project/shmlast/), which uses nucleotides as queries and proteins as target sequences.

The resultings conditional reciprocal best hit pair matrix can be used to obtain pairwise codon alignments, which are further used to calculate synonymous and nonsynonymous substitutions using parallelization. The ka/ks (dN/dS) values can be obtained either via the codon model of [Li WH. (1999)](https://www.ncbi.nlm.nih.gov/pubmed/8433381) implemented in the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) or the model [Yang Z and Nielson R. (2000)](https://www.ncbi.nlm.nih.gov/pubmed/10666704) implemented in [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download).

# Installation

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

### install packages from [Bioconductor](https://www.bioconductor.org/)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### install [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

```
library(devtools)
install_gitlab("mpievolbio-it/crbhits", host = "https://gitlab.gwdg.de")
#install_github("kullrich/CRBHits", build_vignettes = TRUE, dependencies = FALSE)
```

## External tools installation prerequisites

The source code for both prerequisites (LAST, KaKs_Calculator2.0) are forked within [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits). 

- [LAST](http://last.cbrc.jp/) [http://last.cbrc.jp/last-1060.zip](http://last.cbrc.jp/last-1060.zip)

```
#setwd(paste0(find.package("CRBHits"), "/extdata/"))
#system(paste0("unzip last-1060.zip"))
#setwd(paste0(find.package("CRBHits"), "/extdata/last-1060/"))
## see more installation iformation here, if make fails
## last-install-help: http://last.cbrc.jp/doc/last.html
#system("make")
make.last()
```

- [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)

```
#setwd(paste0(find.package("CRBHits"), "/extdata/"))
#system(paste0("tar -xvf KaKs_Calculator2.0.tar.gz"))
#setwd(paste0(find.package("CRBHits"), "/extdata/KaKs_Calculator2.0/src/"))
#system("make clean")
#system("make")
make.KaKs_Calculator2()
```

If you would like to install the latest version of both tools, you need to download the source code and compile again. 

## Vignettes

These vignettes introduce  [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits)

- [CRBHits Basic Vignette](https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/blob/master/vignettes/CRBHitsBasicVignette.Rmd)
- [dNdS Clustering Vignette](https://gitlab.gwdg.de/mpievolbio-it/crbhits/-/blob/master/vignettes/dNdSClusteringVignette.Rmd)

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
- re-write filter for speedup
- self-blast
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

Aubry, S., Kelly, S., Kümpers, B. M., Smith-Unna, R. D., & Hibberd, J. M. (2014). **Deep evolutionary comparison of gene expression identifies parallel recruitment of trans-factors in two independent origins of C4 photosynthesis.** *PLoS genetics*, **10(6)**. [https://doi.org/10.1371/journal.pgen.1004365](https://doi.org/10.1371/journal.pgen.1004365)

Scott, C. (2017). **shmlast: an improved implementation of conditional reciprocal best hits with LAST and Python.** *Journal of Open Source Software*, 2(9), 142. [https://joss.theoj.org/papers/10.21105/joss.00142](https://joss.theoj.org/papers/10.21105/joss.00142)

Kiełbasa, S. M., Wan, R., Sato, K., Horton, P., & Frith, M. C. (2011). **Adaptive seeds tame genomic sequence comparison.** *Genome research*, 21(3), 487-493. [https://doi.org/10.1101/gr.113985.110](https://doi.org/10.1101/gr.113985.110)

Rost, B. (1999). **Twilight zone of protein sequence alignments.** *Protein engineering*, 12(2), 85-94. [https://doi.org/10.1093/protein/12.2.85](https://doi.org/10.1093/protein/12.2.85)

Charif, D., & Lobry, J. R. (2007). **SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.** *In Structural approaches to sequence evolution* (pp. 207-232). Springer, Berlin, Heidelberg. [https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10](https://link.springer.com/chapter/10.1007/978-3-540-35306-5_10)

Li, W. H. (1993). **Unbiased estimation of the rates of synonymous and nonsynonymous substitution.** *Journal of molecular evolution*, 36(1), 96-99. [https://doi.org/10.1007/bf02407308](https://doi.org/10.1007/bf02407308)

Yang, Z., & Nielsen, R. (2000). **Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models.** *Molecular biology and evolution*, 17(1), 32-43. [https://doi.org/10.1093/oxfordjournals.molbev.a026236](https://doi.org/10.1093/oxfordjournals.molbev.a026236)
