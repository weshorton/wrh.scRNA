% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readVsEachMarkers.R
\name{readVsEachMarkers}
\alias{readVsEachMarkers}
\title{Read VS Each Markers}
\usage{
readVsEachMarkers(
  inDir_v,
  dirs_v,
  grep_v = "2x|3xNR|3xR|4x|d80|Ent|NT|PTX|3x|Unt",
  toRead_v = NULL,
  toReadSpop_v = NULL
)
}
\arguments{
\item{inDir_v}{path to the main data directory that contains dirs_v}

\item{dirs_v}{vector of directories to read data from. Should be batch3_lymphoid, batch3_myeloid, batch3_neoplastic (or batch12)}

\item{grep_v}{character string to extract names. Should be "2x|3xNR|3xR|4x|d80|Ent|NT|PTX"}

\item{toRead_v}{used to grep which files to read in within an spop directory}

\item{toReadSpop_v}{used to grep which spop directories to read in}
}
\value{
list of length(dirs_v), each element is a list whose length is the number of treatments run, and each of those
is a list whose length is the number of treatments compared to base treatment. Each element of the list is the data.table
of differential expression results.
}
\description{
Read in directory of vs each find markers results. For each meta population tested (e.g. lymphoid, myeloid, neoplastic),
there is a sub-directory for each subpopulation. For each of those, there is a FindMarkers excel workbook for each treatment.
Each workbook has a sheet for every treatment it's compared to.
Want to read these in so we have a list of list of lists of all of these results
}
