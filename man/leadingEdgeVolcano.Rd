% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leadingEdgeVolcano.R
\name{leadingEdgeVolcano}
\alias{leadingEdgeVolcano}
\title{Leading Edge Volcano}
\usage{
leadingEdgeVolcano(
  data_dt,
  gsea_res,
  pathway_v,
  treat_v,
  lfc_v = 0.5,
  pval_v = 0.05,
  title_v = NULL
)
}
\arguments{
\item{data_dt}{table of DEG. Must be same one that created gsea_res}

\item{gsea_res}{output of GSEA runFGSEA (must be wrangled...change this)}

\item{pathway_v}{vector of the pathway name to make plot for}
}
\value{
ggplot
}
\description{
Make a volcano plot of Leading Edge genes from GSEA results
}
\details{
For a given pathway in gsea_res, grab the leading edge genes and then plot
a volcano plot using the supplied DEG results.
}
