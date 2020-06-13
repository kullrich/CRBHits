---
title: 'CRBHits: Conditional Reciprocal Best Hits in R'
tags:
  - R
  - reciprocal best hit
  - conditional reciprocal best hit
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
[crb-blast](https://github.com/cboursnell/crb-blast) in [R](https://cran.r-project.org/) 
[@team2013r]. The new R package targets ecology, population and evolutionary biologists 
working in the field of comparative genomics.

The Reciprocal Best Hit (RBH) approach is commonly used in bioinformatics to show that two 
sequences evolved from a common ancestral gene. In other words, RBH tries to find orthologous 
protein sequences within and between species. These orthologous sequences can be further 
analysed to evaluate protein family evolution, infer phylogenetic trees and to annotate 
protein function [@altenhoff2019inferring]. The initial sequence search step is classically 
performed with the Basic Local Alignment Search Tool (blast) [@altschul1990basic] and due to 
evolutionary constraints, in most cases protein coding sequences are compared between two 
species. Downstream analysis use the resulting RBH to cluster sequence pairs and build so-
called orthologous groups like e.g. [OrthoFinder](https://github.com/davidemms/OrthoFinder) 
[@emms2015orthofinder] and other tools.

The CRBH algorithm was introduced by @aubry2014deep and builds upon the traditional 
RBH approach to find additional orthologous sequences between two sets of sequences. 
As described earlier [@aubry2014deep; @scott2017shmlast], CRBH uses the sequence search 
results to fit an expect-value (e-value) cutoff given each RBH to subsequently add sequence pairs 
to the list of bona-fide orthologs given their alignment length.

Unfortunately, as mentioned by @scott2017shmlast, the original 
implementation of CRBH ([crb-blast](https://github.com/cboursnell/crb-blast)) lag improved 
blast-like search algorithm to speed up the analysis. As a consequence,
@scott2017shmlast ported CRBH to python [shmlast](https://github.com/camillescott/shmlast), 
while [shmlast](https://github.com/camillescott/shmlast) cannot deal with IUPAC nucleotide
code so far.

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) constitutes a new R package, which 
build upon previous implementations and ports CRBH into the [R](https://cran.r-project.org/) 
environment, which is popular among biologists. 
[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) improve CRBH by additional implemented 
filter steps [@rost1999twilight] and the possibility to apply custom filters.

# Downstream functionalities

Calculating synonymous (dS) and nonsynonymous substitutions (dN) per orthologous sequence 
pair is a common task for evolutionary biologists, since its ratio dN/dS can be used as an 
indicator of selective pressure acting on a protein [@kryazhimskiy2008population]. However, 
this task is computational more demanding and consist of at least two steps, namely 
codon sequence alignment creation and dN/dS calculation. Further, the codon sequence alignment 
step consist of three subtasks, namely coding nucleotide to protein sequence translation, 
pairwise protein sequence alignment calculation and converting the protein 
sequence alignment back into a codon based alignment.

Downstream of CRBH creation, [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) features 
all above mentioned steps and subtasks. [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) 
has the ability to directly create codon alignments within R with the help of the widely used R 
package [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) 
[@pages2017biostrings] (more than 200k downloads per year since 2014). These codon alignments 
can be subsequently used to calculate synonymous and nonsynonymous substitutions per sequence 
pair and is implemented in a multithreaded fashion either via the R package 
[seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) [@charif2007seqinr] or the 
use of an R external tool 
[KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download) [@wang2010kaks_calculator].

# Implementation

Like [shmlast](https://github.com/camillescott/shmlast), 
[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) benefits from the blast-like sequence 
search software [LAST](http://last.cbrc.jp/)[@kielbasa2011adaptive] and plots the fitted model 
of the CRBH e-value based algorithm. In addition, users can filter the hit pairs prior to CRBH 
fitting for other criteria like query coverage, protein identity and/or 
the twilight zone of protein sequence alignments according to 
@rost1999twilight. The implemented filter uses equation 2 [see @rost1999twilight]:

$$f(x_{\text{hit pair}}) = \begin{cases}
100 \text{ , for } L_{\text{hit pair}} < 11 \\
480 * L^{-0.32 * (1 + e^{\frac{-L}{1000}})} \text{ , for } L_{\text{hit pair}} \leq 450 \\
19.5 \text{ , for } L_{\text{hit pair}} > 450
\end{cases}$$

where $x_{\text{hit pair}}$ is the expected protein identity given the alignment length $L_{\text{hit pair}}$. If the actual protein identity of a hit pair exceeds the expected protein identity ($pident_{\text{hit pair}} \geq f(x_{\text{hit pair}})$), it is retained for subsequent CRBH calculation.

In contrast to previous implementations, [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) only take coding nucleotide sequences (CDS) as the query and target inputs. This is due to the downstream functionality of [CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) to directly calculate codon alignments within R, which rely on CDS. The inputs are translated into protein sequences, aligned globally [@smith1981identification] and converted into codon alignments. 

Functions are completely coded in R and only the external prerequisites 
([LAST](http://last.cbrc.jp/) and 
[KaKs_Calculator2.0](https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download)) 
need to be compiled. Further, users can create their own filters before CRBH 
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

The obtained CRBH can also be used to calculate synonymous (dS/ks) and nonsynonymous (dN/ka) substitutions per hit pair using either the model from @li1993unbiased or from @yang2000estimating.

```r
cds1 <- Biostrings::readDNAStringSet(cds1)
cds2 <- Biostrings::readDNAStringSet(cds2)
#get help ?rbh2kaks
cds1.cds2.kaks.Li <- rbh2kaks(cds1.cds2.crbh$crbh.pairs, cds1, cds2,
                              model = "Li", threads = 4)
cds1.cds2.kaks.YN <- rbh2kaks(cds1.cds2.crbh$crbh.pairs, cds1, cds2,
                              model = "YN", threads = 4)
```

Table: Performance comparison for CRBH and dN/dS calculations (Intel Xeon CPU E5-2620 v3 @ 2.40GHz; 3575 hit pairs).\label{tab:performance}

| Number of Threads | 1 | 2 | 4 | 8 |
| - | - | - | - | - | 
| Runtime of CRBH(shmlast) in sec| 38 (s)| 30 (s) | 28 (s) | 28 (s) |
| Runtime of CRBH(CRBHits) in sec| 32 (s)| 26 (s) | 24 (s) | 22 (s) |
| Runtime of kaks.Li in sec| 357 (s)| 167 (s) | 87 (s) | 49 (s) | 
| Runtime of kaks.YN in sec| 474 (s)| 230 (s) | 121 (s) | 63 (s) |

# Conclusions

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) implements CRBH in [R](https://cran.r-project.org/) (see \autoref{fig:crbh}) and also can be used to calculate codon alignment based nucleotide diversities in a multithreaded fashion (see \autoref{tab:performance}).

# Availability

[CRBHits](https://gitlab.gwdg.de/mpievolbio-it/crbhits) is an open source software made available under the MIT license. It can be installed from its gitlab repository using the [devtools](https://devtools.r-lib.org) package.

```r
devtools::install_gitlab("mpievolbio-it/crbhits", 
 host = "https://gitlab.gwdg.de")", build_vignettes = TRUE)
```

The R package website, which contain a detailed HOWTO to install the prerequisites (mentioned above) and package vignettes are availbale at [https://mpievolbio-it.pages.gwdg.de/crbhits](https://mpievolbio-it.pages.gwdg.de/crbhits).

# References
