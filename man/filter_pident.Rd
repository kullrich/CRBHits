% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_pident.R
\name{filter_pident}
\alias{filter_pident}
\title{filter_pident}
\usage{
filter_pident(rbh, pident = 0, inverse = FALSE)
}
\arguments{
\item{rbh}{BLAST-like tabular matrix [mandatory]}

\item{pident}{percent identity [default: 0.0]}

\item{inverse}{specify if filter should keep the removed values [default: FALSE]}
}
\value{
rbh matrix
}
\description{
This function filters BLAST-like tabular output according to protein identity.
}
\examples{
## load crbh data
data(ath_aly_crbh)
dim(ath_aly_crbh$crbh1)
dim(filter_pident(ath_aly_crbh$crbh1, pident = 75))
}
\author{
Kristian K Ullrich
}
