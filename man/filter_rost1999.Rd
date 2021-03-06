% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_rost1999.R
\name{filter_rost1999}
\alias{filter_rost1999}
\title{filter_rost1999}
\usage{
filter_rost1999(rbh, inverse = FALSE)
}
\arguments{
\item{rbh}{BLAST-like tabular matrix [mandatory]}

\item{inverse}{specify if filter should keep the removed values [default: FALSE]}
}
\value{
rbh matrix
}
\description{
This function filters BLAST-like tabular output according to equation 2 of Rost 1999.
}
\examples{
## load crbh data
data(ath_aly_crbh)
dim(ath_aly_crbh$crbh1)
dim(filter_rost1999(ath_aly_crbh$crbh1))
}
\references{
Rost B. (1999). Twilight zone of protein sequence alignments. \emph{Protein Engineering}, \bold{12(2)}, 85-94.
}
\author{
Kristian K Ullrich
}
