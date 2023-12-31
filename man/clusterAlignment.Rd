% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterAlignment.R
\name{clusterAlignment}
\alias{clusterAlignment}
\alias{clusterAlignment-show}
\alias{clusterAlignment-class}
\alias{clusterAlignment-plot}
\alias{show,clusterAlignment-method}
\alias{plot,clusterAlignment-method}
\alias{plot,clusterAlignment,ANY-method}
\title{Data Structure for a collection of all pairwise alignments of GCMS runs}
\usage{
clusterAlignment(
  pD,
  runs = 1:length(pD@rawdata),
  timedf = NULL,
  usePeaks = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{pD}{a \code{peaksDataset} object.}

\item{runs}{vector of integers giving the samples to calculate set of
pairwise alignments over.}

\item{timedf}{list (length = the number of pairwise alignments) of matrices
giving the expected time differences expected at each pair of peaks used
with \code{usePeaks}=\code{TRUE}, passed to \code{peaksAlignment}}

\item{usePeaks}{logical, \code{TRUE} uses \code{peakdata} list, \code{FALSE}
uses \code{rawdata} list for computing similarity.}

\item{verbose}{logical, whether to print out info.}

\item{...}{other arguments passed to \code{peaksAlignment}}
}
\value{
\code{clusterAlignment} object
}
\description{
Store the raw data and optionally, information regarding signal peaks for a
number of GCMS runs
}
\details{
clusterAlignment computes the set of pairwise alignments.
}
\examples{

require(gcspikelite)

# paths and files
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
eluFiles <- dir(gcmsPath, "ELU", full=TRUE)

# read data, peak detection results
pd <- peaksDataset(cdfFiles[1:2], mz=seq(50,550), rtrange=c(7.5,8.5))
pd <- addAMDISPeaks(pd, eluFiles[1:2])

ca <- clusterAlignment(pd, gap=0.5, D=0.05, df=30, metric=1, type=1)

}
\references{
Mark D Robinson (2008).  Methods for the analysis of gas
chromatography - mass spectrometry data \emph{PhD dissertation} University
of Melbourne.
}
\seealso{
\code{\link{peaksDataset}}, \code{\link{peaksAlignment}}
}
\author{
Mark Robinson, Riccardo Romoli
}
\keyword{classes}
