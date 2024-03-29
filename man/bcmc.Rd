% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcmc.R
\name{bcmc}
\alias{bcmc}
\title{BCMC (Biomarker Categorization in Meta-analysis by Concordance)
for biomarker detection and categorization.}
\usage{
bcmc(data.exp, data.clin, meta.de = TRUE, de.lfc = NULL, de.pval = NULL)
}
\arguments{
\item{data.exp:}{The gene expression data. It should be a list in which
every element represents a study. Within each element (study), the data is
supposed to be matrix where each row is a gene and each column is a sample.}

\item{data.clin:}{The clinical data about whether a sample is case or control.
It should be a list in which every element represents a study. Within each
element (study), the data is supposed to be a data frame with only one column
and each row is a sample (same order as in \code{data.exp}). In each data
frame, code control as 0 and case a 1.}

\item{meta.de:}{A logical value indicating whether the \code{MetaDE} should
be run before the BCMC. The default is TRUE. If FALSE, users need to provide
the log fold change and p-value from DE analysis as arguments \code{de.lfc}
and \code{de.pval}.}

\item{de.lfc:}{If \code{meta.de} is FALSE, this should be a matrix of log fold
change from DE analysis. Each row should be a gene and each column should be
a study. Otherwise, the deafult is TRUE and \code{MetaDE} will be run for
DE analysis before the BCMC.}

\item{de.pval:}{If \code{meta.de} is FALSE, this should be a matrix of p-value
from DE analysis. Each row should be a gene and each column should be
a study. Otherwise, the deafult is TRUE and \code{MetaDE} will be run for
DE analysis before the BCMC.}
}
\value{
A list with following components:
\item{lfc}{The matrix of log fold change from DE analysis. Each row is a gene
and each column is a study.}
\item{pval}{The matrix of p-value from DE analysis. Each row is a gene
and each column is a study.}
\item{pos.wp}{The matrix of predicted up-regulated weight pattern. Each row
is a gene and each column is a study.}
\item{neg.wp}{The matrix of predicted down-regulated weight pattern. Each row
is a gene and each column is a study.}
\item{max.wp}{The matrix of predicted dominant weight pattern. Each row
is a gene and each column is a study. This is found by comparing the two BCMC
statistics for up-regualted and down-regulated.}
\item{Rg}{A matrix of BCMC statistcs, as the \eqn{R_g} in the paper. Each row
is a gene. The first column is for up-regulated pattern, the second column
is for down-regulated pattern, and the third column is the larger one of the
previous two statistics. See the paper for more details.}
}
\description{
This function will first run the \code{MetaDE} to perform DE analysis, and
then implement the meta-analysis method BCMC (Biomarker Categorization
in Meta-analysis by Concordance) for biomarker detection and categorization.
}
\examples{
data("SimulDE")
result_bcmc <- bcmc(data.exp=SimulDE$express, data.clin=SimulDE$clin)
names(result_bcmc)
head(result_bcmc$Rg)
head(result_bcmc$pos.wp)
}
