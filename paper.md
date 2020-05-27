---
title: 'CRBHits: Conditional Reciprocal Best Hits in R'
tags:
  - R
  - reciprocal best hit
  - conitional reciprocal best hit
  - codon alignment
  - dN/dS
  - ka/ks
authors:
  - name: Kristian K Ullrich
    orcid: 0000-0003-4308-9626
    affiliation: "1"
affiliations:
 - name: Max Planck Institute for Evolutionary Biology, Scientific IT group, August Thienemann Str. 2, 24306 Plön
   index: 1
date: 5 May 2020
bibliography: paper.bib

---

# Summary

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is a reimplementation of the 
Conditional Reciprocal Best Hit (CRBH) algorithm 
[crb-blast](https://github.com/cboursnell/crb-blast) in 
[R](https://cran.r-project.org/) [@team2013r]. The algorithm was introduced by 
@aubry2014deep and ported to python 
[shmlast](https://github.com/camillescott/shmlast) by @scott2017shmlast, which 
benefits from the blast-like sequence search software 
LAST [@kielbasa2011adaptive]. As described earlier 
[@aubry2014deep; @scott2017shmlast], CRBH builds upon the classical reciprocal 
best hit (RBH) approach to find additional orthologous sequences between two sets of 
sequences by fitting an expect-value cutoff per alignment length. Due to 
evolutionary constraints in most cases protein coding sequences are used and 
compared between two species, whereas downstream analysis use RBH to cluster 
and build orthologous groups like e.g. 
[OrthoFinder](https://github.com/davidemms/OrthoFinder) [@emms2015orthofinder] 
and other tools.

Unfortunately, as mentioned by @scott2017shmlast, the original implementation 
of CRBH ([crb-blast](https://github.com/cboursnell/crb-blast)) lag improved 
blast-like search algorithm to speed up the analysis, while 
[shmlast](https://github.com/camillescott/shmlast) cannot deal with IUPAC 
code so far.

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) constitutes a new R package to 
build upon previous implementations and to improve CRBH by additional filter steps [@rost1999twilight] and the ability to directly create codon alignments 
within R with the help of the R package 
[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) [@pages2017biostrings]. The resulting codon alignments can be subsequently used to calculate synonymous and nonsynonymous substitutions per sequence pair with the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) [@charif2007seqinr] or [KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) [@wang2010kaks_calculator].

# Implementation

Like [shmlast](https://github.com/camillescott/shmlast), 
[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) plots the fitted model of the 
CRBH evalue based algorithm. In addition, users can filter the hit pairs prior 
CRBH fitting for other criteria like query coverage, protein identity and/or 
the twilight zone of protein sequence alignments according to 
@rost1999twilight. The implemented filter uses equation 2 [see @rost1999twilight]:

$$f(x_{\text{hit pair}}) = \begin{cases}
100 \text{ , for } L_{\text{hit pair}} < 11 \\
480 * L^{-0.32 * (1 + e^{\frac{-L}{1000}})} \text{ , for } L_{\text{hit pair}} \leq 450 \\
19.5 \text{ , for } L_{\text{hit pair}} > 450
\end{cases}$$

, where $x_{\text{hit pair}}$ is the expected protein identity given the alignment length $L_{\text{hit pair}}$. If the actual $pident_{\text{hit pair}} \geq f(x_{\text{hit pair}})$, the given hit pair is retained for subsequent CRBH calculation.

In contrast to previous implementations, [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) only take coding nucleotide sequences (CDS) as the query and target inputs and translates them into protein sequences. This is due to the downstream functionality of [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) to directly calculate codon alignments within R, which rely on CDS. Its functions are completely coded in R and only the external prerequisites 
([LAST](http://last.cbrc.jp/) and 
[KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)) 
needs to be compiled. Further, users can create their own filters before CRBH 
calculation.

# Functions and Examples

The following example shows how to obtain CRBH between the coding sequences of *Schizosaccharomyces pombe* [@wood2012pombase] and *Nematostella vectensis* [@apweiler2004protein] by using two URLs as input strings and multiple threads for calculation.

```r
library(CRBHits)
cds1 <- paste0("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/",
               "feature_sequences/cds.fa.gz")
cds2 <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/",
               "Eukaryota/UP000001593_45351_DNA.fasta.gz")
#get help ?cdsfile2rbh
cds1.cds2.crbh <- cdsfile2rbh(cds1, cds2, plotCurve = TRUE, threads = 4)
```

![Accepted secondary reciprocal best hits based on CRBH fitting.\label{fig:crbh}](figure1.png)

The obtained CRBH can further be used to calculate synonymous (dS/ks) and nonsynonymous (dN/ka) substitutions per hit pair using either the model from @li1993unbiased or from @yang2000estimating.

```r
cds1 <- Biostrings::readDNAStringSet(cds1)
cds2 <- Biostrings::readDNAStringSet(cds2)
#get help ?rbh2kaks
cds1.cds2.kaks.Li <- rbh2kaks(cds1.cds2.crbh$crbh.pairs, cds1, cds2,
                              model = "Li", threads = 4)
cds1.cds2.kaks.YN <- rbh2kaks(cds1.cds2.crbh$crbh.pairs, cds1, cds2,
                              model = "YN", threads = 4)
```

Table: Performance comparison for CRBH and dNdS calculations (Intel Xeon CPU E5-2620 v3 @ 2.40GHz; 3575 hit pairs).\label{tab:performance}

| Number of Threads | 1 | 2 | 4 | 8 |
| - | - | - | - | - | 
| Runtime of CRBH in sec| 36 (s)| 29 (s) | 26 (s) | 25 (s) |
| Runtime of kaks.Li in sec| 357 (s)| 167 (s) | 87 (s) | 49 (s) | 
| Runtime of kaks.YN in sec| 474 (s)| 230 (s) | 121 (s) | 63 (s) |

# Conclusions

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) implements CRBH in [R](https://cran.r-project.org/) (see \autoref{fig:crbh}) and further can be used to calculate codon alignment based nucleotide diversities in a multithreaded fashion (see \autoref{tab:performance}).

# Availability

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is an open source software made available under the MIT license. It can be installed from its gitlab repository using the [devtools](https://devtools.r-lib.org) package.

```r
devtools::install_gitlab("mpievolbio-it/crbhits", 
 host = "https://gitlab.gwdg.de")", build_vignettes = TRUE)
```

The R package website, which contain a detailed HOWTO to install the prerequisites (mentioned above) and package vignettes are availbale at [https://mpievolbio-it.pages.gwdg.de/crbhits](https://mpievolbio-it.pages.gwdg.de/crbhits).

# References
