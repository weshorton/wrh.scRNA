% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myClusterSweep.R
\name{myClusterSweep}
\alias{myClusterSweep}
\title{Cluster Parameter Sweep}
\usage{
myClusterSweep(
  seurat_obj,
  embedding_v = seurat_obj@reductions$pca@cell.embeddings,
  reduction_v = "pca",
  reductionName_v = "umap",
  ndims_v = 10,
  res_v = seq(from = 0.1, to = 1, by = 0.1),
  verbose_v = F,
  indPlots_v = F,
  plotDir_v = NULL,
  name_v = "",
  print_v = T
)
}
\arguments{
\item{embedding_v}{embeddings for chosen dimensional reduction. Usually: seruat_obj@reductions$pca@cell.embeddings}

\item{reduction_v}{which reduction to use. Default is pca (which goes with default embedding). Be sure to change both!}

\item{reductionName_v}{name to call reduction on seurat object. Default is "umap", other values likely "umapR" or "umapH" for rliger/harmony}

\item{ndims_v}{number of dimensions to use. Default is 10.}

\item{res_v}{vector of resolution values to pass to FindClusters()}

\item{verbose_v}{logical indicating whether to print seurat function messages or not.}

\item{indPlots_v}{logical indicating whether to print and/or save the individual QC plots}

\item{plotDir_v}{path to output directory for plots. If NULL, will print to console.}

\item{name_v}{Use as a prefix for output plots. Sample1 for example.}

\item{print_v}{logical. Should main plots be printed, even if they're also saved? Does not change individual plots}

\item{seruat_obj}{A seurat object with dimensional reduction embeddings calculated (PCA)}
}
\value{
Either save plots to file or print to console.
}
\description{
Run a sweep on different resolutions for clustering. Output plots to assess quality of clustering.
Adapted from code provided by Nick Calistri.
}
