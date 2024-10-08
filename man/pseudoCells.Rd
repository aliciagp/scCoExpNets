% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scGCNs_creation.R
\name{pseudoCells}
\alias{pseudoCells}
\title{pseudoCells - It creates a pseudo-cells matrix and the corresponding pseudo-cells metadata table.}
\usage{
pseudoCells(
  exprData,
  covs,
  save = F,
  outputDir = paste0(getwd(), "/pseudoCells/")
)
}
\arguments{
\item{exprData}{it represents the gene expression matrix of a specific cell type cluster}

\item{covs}{it represents the metadata table of a specific cell type cluster}

\item{save}{if we want to save the created files}

\item{outputDir}{if save=TRUE, in which directory we would like to save the files created}
}
\value{
it returns a list where the first element represents the pseudo-cells matrix and the second one represents the pseudo-cells metadata table
}
\description{
pseudoCells - It creates a pseudo-cells matrix and the corresponding pseudo-cells metadata table.
}
