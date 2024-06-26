% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkForGene.R
\name{checkForGene}
\alias{checkForGene}
\title{Check For Gene}
\usage{
checkForGene(deg_lslslslsdt, gene_v, pval_v = NULL)
}
\arguments{
\item{deg_lslslslsdt}{either a vsEach or a global DEG output.}

\item{gene_v}{gene to check}

\item{pval_v}{optional numeric value to subset the deg tables for significant hits of gene_v only}
}
\value{
data.table with one row per hit of gene_v in deg_lslslslsdt
}
\description{
Check all of the deg tables in either the vsEach or global outputs for a given gene
}
