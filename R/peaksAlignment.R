##' Compression method for peaksAlignment object
##'
##' Compression method for peaksAlignment object
##' @param object peaksAlignment
##' @param verbose logical
##' @param ... 
##' @author MR
##' @import methods
##' @importFrom methods setMethod new
##' @keywords internal
setMethod("compress","peaksAlignment",
          function(object,verbose = TRUE, ...) {
              if(object@compressed)
              {
                  if (verbose){
                      cat("[compress.peaksAlignment] Already compressed.\n")}
                  return(object)
              }
              #if (verbose)
              #  before<-sum(ll(object)$KB)
              object@r <- as.matrix.csc(1-object@r)
              #if (verbose) {
              #  after<-sum(ll(object)$KB)
              #  cat("[compress.peaksAlignment] Object is",(before-after),"KB
              #  smaller.\n") 
              #}
              object@compressed <- TRUE
              new("peaksAlignment", object)
          })

##' Decompression method for peaksAlignment object
##'
##' Decompression method for peaksAlignment object
##' @param object peaksAlignment object
##' @param verbose dummy
##' @param ... dummy
##' @author MR
##' @import methods
##' @importFrom methods setMethod new
##' @keywords internal
setMethod("decompress", "peaksAlignment",
          function(object, verbose = TRUE, ...) {
              if(!(object@compressed)) {
                  if(verbose)
                  {
                      cat("[decompress.peaksAlignment] Already decompressed.
\n")
                  }
                  return(object)
              }
              ## if (verbose)
              ##  {before<-sum(ll(object)$KB)}
              object@r <- 1-(as.matrix(object@r))
              ## if (verbose) {
              ##   after<-sum(ll(object)$KB)
              ##   cat("[decompress.peaksAlignment] Object
              ##   is",(after-before),"KB larger.\n") 
              ## }
              object@compressed <- FALSE
              new("peaksAlignment", object)
          })



#' Data Structure for pairwise alignment of 2 GCMS samples
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' peaksAlignment is a hold-all data structure of the raw and peak detection
#' data. 
#' @param object peaksAlignment object
#' @author Mark Robinson, Riccardo Romoli
#' @export
#' @noRd
setMethod("show","peaksAlignment",
          function(object){
              cat("An object of class \"", class(object), "\"\n", sep="")
              dm <- dim(object@r)
              cat(dm[1], "peaks against", dm[2],"peaks:",
                  nrow(object@v$match), "matched\n")
          })



#' Data Structure for pairwise alignment of 2 GCMS samples
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' peaksAlignment is a hold-all data structure of the raw and peak detection
#' data.
#'
#' @name peaksAlignment-class
#' @aliases peaksAlignment-class peaksAlignment-show peaksAlignment-plot
#' peaksAlignment show,peaksAlignment-method plot,peaksAlignment-method
#' plot,peaksAlignment,ANY-method
#' @param d1 matrix of MS intensities for 1st sample (if doing a peak
#' alignment, this contains peak apexes/areas; if doing a profile alignment,
#' this contains scan intensities.  Rows are m/z bins, columns are
#' peaks/scans.
#' @param d2 matrix of MS intensities for 2nd sample
#' @param t1 vector of retention times for 1st sample
#' @param t2 vector of retention times for 2nd sample
#' @param gap gap penalty for dynamic programming algorithm. Not used if
#' \code{type=2}
#' @param D time window (on same scale as retention time differences,
#' \code{t1}
#' and \code{t2}. Default scale is seconds.)
#' @param timedf list (length = the number of pairwise alignments) of matrices
#' giving the expected time differences expected at each pair of peaks used
#' with \code{usePeaks}=\code{TRUE}.
#' @param df integer, how far from the diagonal to go to calculate the
#' similarity of peaks. Smaller value should run faster, but be careful not to
#' choose too low.
#' @param verbose logical, whether to print out info.
#' @param usePeaks logical, \code{TRUE} uses \code{peakdata} list,
#' \code{FALSE}
#' uses \code{rawdata} list for computing similarity.
#' @param compress logical, whether to compress the similarity matrix into a
#' sparse format.
#' @param metric numeric, different algorithm to calculate the similarity
#' matrix between two mass spectrum. \code{metric=1} call
#' \code{normDotProduct()}; \code{metric=2} call \code{ndpRT()};
#' \code{metric=3} call \code{corPrt()}
#' @param type numeric, two different type of alignment function
#' @param penality penalization applied to the matching between two mass
#' spectra if \code{(t1-t2)>D}
#' @return \code{peaksAlignment} object
#' @author Mark Robinson, Riccardo Romoli
#' @seealso \code{\link{peaksDataset}}, \code{\link{clusterAlignment}}
#' @references Mark D Robinson (2008).  Methods for the analysis of gas
#' chromatography - mass spectrometry data \emph{PhD dissertation} University
#' of Melbourne.
#' @keywords classes
#' @examples
#' 
#' ## see clusterAlignment, it calls peaksAlignment
#' 
#' ## Not Run:
#' files <- list.files(path = paste(find.package("gcspikelite"), "data",
#'                     sep = "/"),"CDF", full = TRUE)
#' data <- peaksDataset(files[1:2], mz = seq(50, 550), rtrange = c(7.5, 8.5))
#' ## create settings object
#' mfp <- xcms::MatchedFilterParam(fwhm = 10, snthresh = 5)
#' cwt <- xcms::CentWaveParam(snthresh = 3, ppm = 3000, peakwidth = c(3, 40),
#'  prefilter = c(3, 100), fitgauss = FALSE, integrate = 2, noise = 0,
#'  extendLengthMSW = TRUE, mzCenterFun = "wMean")
#' data <- addXCMSPeaks(files[1:2], data, settings = mfp, minintens = 100,
#'  multipleMatchedFilter = FALSE, multipleMatchedFilterParam =
#'  list(fwhm = c(5, 10, 20), rt_abs = 3, mz_abs = 0.1))
#' data
#' plotChrom(data, rtrange=c(7.5, 10.5), runs=c(1:2))
#' 
#' ## align two chromatogram
#' pA <- peaksAlignment(data@peaksdata[[1]], data@peaksdata[[2]],
#'                      data@peaksrt[[1]], data@peaksrt[[2]], D = 50,
#'                      metric = 3, compress = FALSE, type = 2, penality = 0.2)
#' 
#' plotAlignment(pA)
#' pA@v$match
#' 
#' par(mfrow=c(2,1))
#' plot(data@peaksdata[[1]][,15], type = 'h', main = paste(data@peaksrt[[1]][[15]]))
#' plot(data@peaksdata[[2]][,17], type = 'h',
#'      main = paste(data@peaksrt[[2]][[17]]))
#' ## End (Not Run)
#'
#' @export
peaksAlignment <- function(d1, d2, t1, t2, gap = 0.5, D = 50,
                           timedf = NULL, df = 30, verbose = TRUE, 
                           usePeaks = TRUE, compress = TRUE, metric = 2,
                           type = 2, penality = 0.2) {

    ## r <- switch(metric,
    ##             normDotProduct(d1,d2,t1,t2,D=D,
    ##                            df=df+abs(ncol(d1)-ncol(d2)),
    ##                            timedf=timedf,verbose=verbose), # metric=1 
    ##             cos.ndp.rank(d1,d2,t1,t2,D=D,timedf=timedf,
    ##                          verbose=verbose), # metricic=2
    ##             euclidean(d1,d2,t1,t2,D=D,timedf=timedf,
    ##                       verbose=verbose),# metric=3
    ##             soai(d1,d2,t1,t2,D=D,timedf=timedf,
    ##                  verbose=verbose), # metric=4
    ##             manhattan(d1,d2,t1,t2,D=D,timedf=timedf,
    ##                       verbose=verbose), # metric=5
    ##             window.metric.c(d1,d2,mzind=mzind,df=df,w=w) # metric=6
    ##             )
                                                         

    ## if(metric == 1) # why?
    ## {
    ##     D <- D/100 
    ## }
    r <- switch(metric,
                normDotProduct(d1, d2, t1, t2, D = D,
                               df = df + abs(ncol(d1) - ncol(d2)),
                               timedf = timedf, verbose = verbose),
                ndpRT(s1 = d1, s2 = d2, t1, t2, D = D),
                corPrt(d1, d2, t1, t2, D = D, penality = penality)
                )
    
    r[is.nan(r)] <- 1 ## remove NaN

    if (type == 1) {
        if(verbose) 
            cat("[peaksAlignment] Comparing", ncol(d1), "peaks to", 
                ncol(d2), "peaks -- gap=", gap, "D=", D, ', metric=',
                metric, ', type=', type, "\n")
        v <- dp(r, gap = gap, verbose = verbose) # dynamic programming
    }
    
    if (type == 2) {
        if(verbose) 
            cat("[peaksAlignment] Comparing", ncol(d1), "peaks to", 
                ncol(d2), "peaks -- D=", D, "seconds,", 'metric=',
                metric, ', type=', type, '\n')
        v <- dynRT(S = r) # RR, modified to be used with normDotProduct()
    }
    
    v$match <- v$match[!is.na(v$match[,2]), ] # remove non-matched peaks
    
    sim <- 0
    for(i in 1:nrow(v$match)) {
        sim <- sim + r[v$match[i, 1], v$match[i, 2]]#
    }
    sim <- sim / nrow(v$match)
    if(verbose) {
        cat("[peaksAlignment] ", nrow(v$match), "matched.  Similarity=",
        sim, "\n")
        }
    object <- new("peaksAlignment", v = v, r = r, dist = sim, 
                  compressed = FALSE, gap = gap, D = D)
    if(compress) {
        compress(object)
    } else {
        object
    }
}


#' Plotting functions for GCMS data objects
#'
#' Plot an object of \code{\linkS4class{peaksAlignment}}
#'
#' The similarity matrix is plotted and optionally, the set of matching peaks.
#' \code{clusterAlignment} objects are just a collection of all pairwise
#' \code{peakAlignment} objects.
#'
#' @title plotAlignment
#' @param object a \code{clusterAlignment} object
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param plotMatches logical, whether to plot matches
#' @param matchPch match plotting character
#' @param matchLwd match line width
#' @param matchCex match character expansion factor
#' @param matchCol match colour
#' @param col vector of colours for colourscale
#' @param breaks vector of breaks for colourscale
#' @param ... further arguments passed to \code{image}
#' @return plot an object of class \code{\linkS4class{peaksAlignment}}
#' @author Mark Robinson
#' @seealso  \code{\link{peaksAlignment}} \code{\link{plotAlignment}}
#' @references Mark D Robinson (2008).  Methods for the analysis of gas
#' chromatography - mass spectrometry data \emph{PhD dissertation} University
#' of Melbourne.
#' @keywords classes
#' @examples
#'
#' require(gcspikelite)
#' files <- list.files(path = paste(find.package("gcspikelite"), "data",
#'                     sep = "/"),"CDF", full = TRUE)
#' data <- peaksDataset(files[1:2], mz = seq(50, 550), rtrange = c(7.5, 8.5))
#' ## create settings object
#' mfp <- xcms::MatchedFilterParam(fwhm = 10, snthresh = 5)
#' cwt <- xcms::CentWaveParam(snthresh = 3, ppm = 3000, peakwidth = c(3, 40),
#'  prefilter = c(3, 100), fitgauss = FALSE, integrate = 2, noise = 0,
#'  extendLengthMSW = TRUE, mzCenterFun = "wMean")
#' data <- addXCMSPeaks(files[1:2], data, settings = mfp, minintens = 100,
#'  multipleMatchedFilter = FALSE, multipleMatchedFilterParam =
#'  list(fwhm = c(5, 10, 20), rt_abs = 3, mz_abs = 0.1))
#' data
#' ## image plot
#' plotChrom(data, rtrange = c(7.5,8.5), plotPeaks = TRUE, plotPeakLabels =TRUE)
#'
#' ## align two chromatogram
#' pA <- peaksAlignment(data@peaksdata[[1]], data@peaksdata[[2]],
#'                      data@peaksrt[[1]], data@peaksrt[[2]], D = 50,
#'                      compress = FALSE, type = 1, metric = 1,
#'                      gap = 0.5)
#' plotAlignment(pA)
#'
#' @importFrom graphics image points
#' @export
setMethod("plotAlignment", "peaksAlignment",
          function(object, xlab = "Peaks - run 1",
                   ylab = "Peaks - run 2", plotMatches = TRUE,
                   matchPch = 19, matchLwd = 3, matchCex = 0.5,
                   matchCol = "black",
                   ## col=colorpanel(50,"black", "blue","white"),# 
                   col = colorpanel(50,'white', "green","navyblue"),
                   breaks = seq(0, 1, length = 51), ...)
          {
              if(object@compressed)
              {
                  object <- decompress(object, ...)
              }
              r <- object@r
              image(1:nrow(r), 1:ncol(r), r, col = col, xlab = xlab,
                    ylab = ylab, ...)
              if(plotMatches)
              {
        points(object@v$match, pch = matchPch, col = matchCol,
               cex = matchCex, type = "b", lwd = matchLwd)
              }
          }
        )

