% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/globalDeletion.R
\name{globalDeletion}
\alias{globalDeletion}
\title{globalDeletion}
\usage{
globalDeletion(dna)
}
\arguments{
\item{dna}{\code{DNAStringSet}}
}
\description{
This function returns a \code{DNAStringSet} reduced by all
sites containing any gaps ("-", "+", ".") or missing ("N") sites.
}
\author{
Kristian K Ullrich
}
