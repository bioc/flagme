setMethod("compress","clusterAlignment",
  function(object, verbose = TRUE, ...){
    for(i in 1:length(object@alignments)){
      pA <- object@alignments[[i]]
      object@alignments[[i]] <- compress(pA, verbose = FALSE)
    }
  new("clusterAlignment", object)
})

setMethod("decompress","clusterAlignment",
  function(object, verbose = TRUE, ...){
    for(i in 1:length(object@alignments)){
      pA <- object@alignments[[i]]
      object@alignments[[i]] <- decompress(pA, verbose = FALSE)
    }
  new("clusterAlignment", object)
})





#' Data Structure for a collection of all pairwise alignments of GCMS runs
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' clusterAlignment computes the set of pairwise alignments.
#' 
#' @aliases clusterAlignment clusterAlignment-show clusterAlignment-class
#' clusterAlignment-plot show,clusterAlignment-method
#' plot,clusterAlignment-method plot,clusterAlignment,ANY-method
#' @param pD a \code{peaksDataset} object.
#' @param runs vector of integers giving the samples to calculate set of
#' pairwise alignments over.
#' @param timedf list (length = the number of pairwise alignments) of matrices
#' giving the expected time differences expected at each pair of peaks used
#' with \code{usePeaks}=\code{TRUE}, passed to \code{peaksAlignment}
#' @param usePeaks logical, \code{TRUE} uses \code{peakdata} list, \code{FALSE}
#' uses \code{rawdata} list for computing similarity.
#' @param verbose logical, whether to print out info.
#' @param ... other arguments passed to \code{peaksAlignment}
#' @return \code{clusterAlignment} object
#' @author Mark Robinson, Riccardo Romoli
#' @seealso \code{\link{peaksDataset}}, \code{\link{peaksAlignment}}
#' @references Mark D Robinson (2008).  Methods for the analysis of gas
#' chromatography - mass spectrometry data \emph{PhD dissertation} University
#' of Melbourne.
#' @keywords classes
#' @examples
#' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
#' cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
#' eluFiles <- dir(gcmsPath, "ELU", full=TRUE)
#' 
#' # read data, peak detection results
#' pd <- peaksDataset(cdfFiles[1:2], mz=seq(50,550), rtrange=c(7.5,8.5))
#' pd <- addAMDISPeaks(pd, eluFiles[1:2])
#' 
#' ca <- clusterAlignment(pd, gap=0.5, D=0.05, df=30, metric=1, type=1)
#'
#' @importFrom stats as.dist hclust
#' @export clusterAlignment
clusterAlignment <- function(pD, runs = 1:length(pD@rawdata),
                             timedf = NULL, usePeaks = TRUE, 
                             verbose = TRUE, ...) {
    n <- length(runs)
    if (usePeaks)
        nr <- length(pD@peaksdata)
    else
        nr <- length(pD@rawdata)
    alignments <- vector("list", n * (n - 1) / 2)
    aligned <- matrix(-1, nr, nr)
    colnames(aligned) <- names(pD@rawdata)
    rownames(aligned) <- names(pD@rawdata)
    dist <- matrix(0, n, n)
    colnames(dist) <- names(pD@rawdata)[runs]
    rownames(dist) <- names(pD@rawdata)[runs]
    count <- 0
    for (i in 1:( n - 1)) {
        run.i <- runs[i]
        for (j in (i + 1):n) {
            run.j <- runs[j]
            count <- count + 1
            if(verbose) {
                cat("[clusterAlignment] Aligning",
                    names(pD@rawdata)[run.i], "to",
                    names(pD@rawdata)[run.j], "\n")
            }
            if (usePeaks) {
              alignments[[count]] <-
                peaksAlignment(pD@peaksdata[[run.i]],
                               pD@peaksdata[[run.j]],
                               pD@peaksrt[[run.i]],
                               pD@peaksrt[[run.j]],
                               usePeaks = usePeaks,
                               timedf = timedf[[count]],
                               verbose = verbose, ...)
            } else {
                alignments[[count]] <-
                    peaksAlignment(pD@rawdata[[run.i]],
                                   pD@rawdata[[run.j]],
                                   pD@rawrt[[run.i]],
                                   pD@rawrt[[run.j]],
                                   usePeaks = usePeaks, timedf = NULL,
                                   verbose = verbose, ...)
            }
            aligned[runs[i], runs[j]] <- aligned[runs[j], runs[i]] <- count
            dist[j, i] <- dist[i, j] <- alignments[[count]]@dist
          }
  }
    merge <- hclust(as.dist(dist), method = "average")$merge
    merge.copy <- merge
    for (i in 1:length(runs)) {
      merge[which(merge.copy == (-i))] <- (-runs[i])
      }
    new("clusterAlignment", runs = runs, aligned = aligned,
        gap = alignments[[1]]@gap, D = alignments[[1]]@D, dist = dist,
        alignments = alignments, merge = merge)
}




#' Data Structure for a collection of all pairwise alignments of GCMS runs
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' clusterAlignment computes the set of pairwise alignments.
##' @param object a clusterAlignment object
##' @author Mark Robinson, Riccardo Romoli
##' @noRd
setMethod("show","clusterAlignment",
function(object) {
   cat("An object of class \"",class(object),"\"\n",sep="")
   cat("Pairwise distance matrix\n")
   print(object@dist)
   #cat(length(object@mz), "m/z bins - range: (",range(object@mz),")\n",sep=" ")
   #cat("scans:",sapply(object@rawdata,ncol),"\n",sep=" ")
   #cat("peaks:",sapply(object@peaksdata,ncol),"\n",sep=" ")
})



#' Plotting functions for GCMS data objects
#'
#' For \code{clusterAlignment} objects, the similarity matrix is plotted and
#' optionally, the set of matching peaks.  \code{clusterAlignment} objects are
#' just a collection of all pairwise \code{peakAlignment} objects.
##' @title plotClustAlignment
##' @param object \code{clusterAlignment} object.
##' @param alignment the set of alignments to plot
##' @param ... further arguments passed to \code{image}. See
##' also \code{plotAlignment}
##' @return plot the pairwise alignment
##' @author Mark Robinson
##' @seealso \code{\link{plotAlignment}}
##' @references Mark D Robinson (2008).  Methods for the analysis of gas
##' chromatography - mass spectrometry data \emph{PhD dissertation} University
##' of Melbourne.
##' @keywords classes
##' @examples
##' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
#' cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
#' eluFiles <- dir(gcmsPath, "ELU", full=TRUE)
#' 
#' # read data, peak detection results
#' pd <- peaksDataset(cdfFiles[1:2], mz=seq(50,550), rtrange=c(7.5,8.5))
#' pd <- addAMDISPeaks(pd, eluFiles[1:2])
#' 
#' ca <- clusterAlignment(pd, gap=0.5, D=0.05, df=30, metric=1, type=1)
##' plotClustAlignment(ca, run = 1)
##' plotClustAlignment(ca, run = 2)
##' plotClustAlignment(ca, run = 3)
##' 
##' @importFrom graphics plot
##' @export
setMethod("plotClustAlignment", "clusterAlignment",
          function(object, alignment = 1, ...) {
              rn <- rownames(object@aligned)
              for (i in alignment) {
                  ind <- which(object@aligned == i, arr.ind = TRUE)[2, ]
                  plot(## object@alignments[[i]],
                      object@alignments[[i]]@v$match,
                      main = paste("D=", object@D, " gap=", object@gap,
                                   sep = ""),
                      xlab = paste("Peaks ", rn[ind[1]], sep = " - "),
                      ylab = paste("Peaks ", rn[ind[2]], sep = " - "),
                      ...)
              }
          }
          )
