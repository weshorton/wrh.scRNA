% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myLISI.R
\name{myLISI}
\alias{myLISI}
\title{Calc and Plot LISI}
\usage{
myLISI(seurat_obj, reduction_v, metaCols_v, name_v)
}
\arguments{
\item{seurat_obj}{seurat object to use}

\item{reduction_v}{name of reduction to use for embeddings}

\item{metaCols_v}{vector of column names to extract from metadata. First must be the grouping variable
(e.g. batchID); second must be the cell classifier variable}

\item{name_v}{name to prepend to plots}
}
\value{
list of ggplot objects
}
\description{
Wrapper function to run compute_lisi and make plots
}
\details{
This still needs some work. Right now must have both meta columns (batch and cell class)
}
