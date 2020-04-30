---
title: "dNdS Clustering Vignette"
author: "Kristian K Ullrich"
date: "2020-04-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{dNdS Clustering Vignette}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---



1. [Conditional reciprocal best hit pairs](#crbhpairs)

To calculate conditional reciprocal best hit pairs between two species, one can directly use an URL to access the coding sequences and calculate the conditional reciprocal best hit pair matrix (crbh).

As an example, here the coding sequences from *Arabidopsis thaliana* and *Arabidopsis lyrata* are used as input sequences from the FTP server from [EnsemblPlants](https://plants.ensembl.org/index.html).

The [release-47](ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/) is used.


```r
library(CRBHits)
athfile <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-47/
fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz"
alyfile <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-47/
fasta/arabidopsis_lyrata/cds/Arabidopsis_lyrata.v.1.0.cds.all.fa.gz"
```

The blast-like software [LAST](http://last.cbrc.jp/) is used to compare the translated coding sequences against each other and output a blast-like output table including the query and target length. Taking the length of the obtained pairwise protein alignment one can calculate for each hit pair the query coverage as $\frac{alignment length}{query length}$.

Next toquery coverage the hit pairs will be additionally filtered for the twilight zone of protein sequence alignments according to [Rost B. (1999)](https://academic.oup.com/peds/article/12/2/85/1550637).

The implemented filter uses `equation2` of the mentioned paper

\begin{equation}
u(x) = 
\begin{cases} 
\exp{(x)} & \text{if } x \geq 0 \\
1         & \text{if } x < 0
\end{cases}
(\#eq:piece)
\end{equation}

$$f(k) = \begin{cases}{n \choose k} p^{k} (1-p)^{n-k};\end{cases}$$.



2. [Reciprocal best hit matrix](#rbhmatrix)
3. [Conditional reciprocal best hit matrix](#crbhmatrix)
