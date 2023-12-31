% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePairs.R
\name{makePairs}
\alias{makePairs}
\title{Make Pairs}
\usage{
makePairs(vector_v, pairs_lsv = NULL, counter_v = 1)
}
\arguments{
\item{vector_v}{vector of elements to combine pairwise}

\item{counter_v}{counter for recursion}

\item{paris_lsv}{list of pairs fed back in recursively}
}
\value{
list of pairs. Each element of vector_v will be a list element whose values are the other elements of vector_v.
}
\description{
Get all unique pairs of values from a vector.
(order doesn't matter)
}
