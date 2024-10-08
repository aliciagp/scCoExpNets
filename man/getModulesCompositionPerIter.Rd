% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scGCNs_visualization.R
\name{getModulesCompositionPerIter}
\alias{getModulesCompositionPerIter}
\title{getModulesCompositionPerIter  - It returns the genes that made up each module in each iteration.
It is an auxiliary function to create a sunburst plot, sankey plot or enrichment plot.}
\usage{
getModulesCompositionPerIter(nets.dir)
}
\arguments{
\item{nets_dir}{the path where we can find all the networks belonging to the same cell type}
}
\value{
a data frame containing the number of genes that remain together at the end after successive pseudo-cells iterations
using all possible combinations of modules from different iterations
}
\description{
getModulesCompositionPerIter  - It returns the genes that made up each module in each iteration.
It is an auxiliary function to create a sunburst plot, sankey plot or enrichment plot.
}
