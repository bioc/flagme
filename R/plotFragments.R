#' plotFrags
#'
#' Plot the mass spectra from the profile matrix
#'
#' Plot the deconvoluted mass spectra from the profile matrix
#'
#' @param object an object of class "peaksDataset" where to keep the mass
#' spectra; both abundance (y) than m/z (x)
#' @param sample character, the sample from were to plot the mass spectra
#' @param specID numerical, a vector containing the index of the spectra to be
#' plotted.
#' @param normalize logical, if TRUE normalize the intensity of the mass peak
#' to 100, the most abundant is 100% and the other peaks are scaled
#' consequetially
#' @param ... other parameter passed to the plot() function
#' @author riccardo.romoli@@unifi.it
#' @examples
#'
#' files <- list.files(path = paste(find.package("gcspikelite"), "data",
#'                     sep = "/"),"CDF", full = TRUE)
#' data <- peaksDataset(files[1:2], mz = seq(50, 550), rtrange = c(7.5, 8.5))
#' ## create settings object
#' mfp <- xcms::MatchedFilterParam(fwhm = 10, snthresh = 5)
#' cwt <- xcms::CentWaveParam(snthresh = 3, ppm = 3000, peakwidth = c(3, 40),
#'  prefilter = c(3, 100), fitgauss = FALSE, integrate = 2, noise = 0,
#'  extendLengthMSW = TRUE, mzCenterFun = "wMean")
#' data <- addXCMSPeaks(files[1:2], data, settings = mfp)
#' data
#' ## align two chromatogram
#' pA <- peaksAlignment(data@peaksdata[[1]], data@peaksdata[[2]],
#'                      data@peaksrt[[1]], data@peaksrt[[2]], D = 50,
#'                      metric = 3, compress = FALSE, type = 2, penality = 0.2)
#' pA@v$match
#' ## plot the mass spectra
#' par(mfrow=c(2,1))
#' plotFrags(object=data, sample=1, specID=10)
#' plotFrags(object=data, sample=2, specID=12)
#'
#' @export plotFrags
plotFrags <- function(object, sample, specID, normalize = TRUE, ...) {
    sp <- as.numeric(specID)
    if (length(sp) > 0) {
        y <- object@peaksdata[[sample]][, sp]
        if (normalize) {
            y <- sapply(y, function(x) {
                i <- which.max(y)
                (100 * x) / y[i]}
                )
                }
        plot(y, type = "h", main = paste("Sample",
            names(object@peaksdata[sample]), "RT",
            object@peaksrt[[sample]][sp], "min"), xlab = "m/z",
             ylab =
                 if (normalize) {
                     "Rel. Abundance"
                 } else {
                     "Abs. Abundance"
                 }, ...)
        } else {
        stop(paste("The spectrum is not present in the sample",
        object@files[sample], '\n'))
        }
}






#' plotAlignedFrags
#'
#' Plot the aligned mass spectra
#'
#' Plot the deconvoluted and aligned mass spectra collected using gatherInfo()
#'
#' @param object where to keep the mass range of the experiment
#' @param outList where to keep the mass spectra; both abundance than m/z
#' @param specID a vector containing the index of the spectra to be plotted. Is
#' referred to outList
#' @param fullRange if TRUE uses the mass range of the whole experiment,
#' otherwise uses only the mass range of each plotted spectum
#' @param normalize if TRUE normalize the intensity of the mass peak to 100,
#' the most abundant is 100\% and the other peaks are scaled consequetially
#' @param ... further arguments passed to the ‘plot’ command
#' @author Riccardo Romoli (riccardo.romoli@@unifi.it)
#' @keywords gatherInfo() plot()
#' @examples
#'
#' files <- list.files(path = paste(find.package("gcspikelite"), "data",
#'                     sep = "/"),"CDF", full = TRUE)
#' data <- peaksDataset(files[1:4], mz = seq(50, 550), rtrange = c(7.5, 8.5))
#' ## create settings object
#' mfp <- xcms::MatchedFilterParam(fwhm = 10, snthresh = 5)
#' cwt <- xcms::CentWaveParam(snthresh = 3, ppm = 3000, peakwidth = c(3, 40),
#'  prefilter = c(3, 100), fitgauss = FALSE, integrate = 2, noise = 0,
#'  extendLengthMSW = TRUE, mzCenterFun = "wMean")
#' data <- addXCMSPeaks(files[1:4], data, settings = mfp)
#' data
#' ## multiple alignment
#' ma <- multipleAlignment(data, c(1,1,2,2), wn.gap = 0.5, wn.D = 0.05,
#'  bw.gap = 0.6, bw.D = 0.2, usePeaks = TRUE, filterMin = 1, df = 50,
#'  verbose = TRUE, metric = 2, type = 2)
#'
#' ## gather apex intensities
#' gip <- gatherInfo(data, ma)
#' gip[[33]]
#' plotAlignedFrags(object = data, outList = gip, specID = 33)
#'
#' @export plotAlignedFrags
plotAlignedFrags <- function(object, outList, specID, fullRange = TRUE,
                             normalize = TRUE, ...) {
  specID <- as.numeric(specID)
  mz <- outList[[specID]]$mz
  abundance <- outList[[specID]]$data
  if (normalize) {
    for (i in 1:ncol(abundance)) {
      if (is.na(sum(abundance[,i])))
        next
      abundance[, i] <- 100 * abundance[, i] /
        abundance[which.max(abundance[, i]), i]
    }
  }

## set the plot grid
  specnum <- table(apply(abundance, 2, sum)) # count number of massSpec different from NA
  if (length(specnum) < 0) {
    paste0("The spec is not present!")
  } else if (length(specnum) <= 3 & length(specnum) >= 1) {
    par(mfrow = c(3, 1))
  } else if(length(specnum) <= 6 & length(specnum) > 3) {
    par(mfrow = c(3, 2))
  } else if(length(specnum) <= 9 & length(specnum) > 6) {
    par(mfrow = c(3, 3))
  } else if(length(specnum) <= 12 & length(specnum) > 9) {
    par(mfrow = c(4, 3))
  } else if(length(specnum) <= 16 & length(specnum) > 12) {
    par(mfrow = c(4, 4))
  } else if(length(specnum) > 16) {
    paste0("Too many spectra to plot. I will plot a random selection of the
spec", specID, " from among the sample list")
    par(mfrow = c(4, 4))
    idx <- sample(1:ncol(abundance), size = 16)
    abundance <- abundance[,idx]
  }
## do the plots
for (i in 1:ncol(abundance)){
  if (is.na(apply(abundance, 2, sum)[i]))
    next
  plot(mz, abundance[, i], type = 'h', main = paste('Retention Time',
                    round(mean(outList[[specID]]$rt[i], na.rm = TRUE),
                          digits = 2), 'min - Sample ', basename(object@files[i])),
       xlab = 'm/z', ylab = if (normalize) {
                              'Rel. Abundance'
                            } else {
                              'Abs. Abundance'
                            },
       xlim = if(fullRange) {
                range(object@mz)
              }
       )
}
}



## original by MR
##
## .plotFragments<-function(peaks,pk=1,cols=NULL,ltys=NULL,TRANSFUN=log2,...) {
##   if(is.null(cols))
##     cols<-rep(c("black","green","blue","red"),each=ncol(peaks[[pk]]$data)/4)
##   if(is.null(ltys))
##      ltys<-rep(1,length(peaks))
##   pr<-paste(round(range(peaks[[pk]]$rt,na.rm=TRUE),2),collapse="-")
##   matplot(peaks[[pk]]$mz,TRANSFUN(peaks[[pk]]$data),type="l",col=cols,lty=ltys,xlab="m/z",ylab="Intensity",main=paste("index:",pk,"rt:",pr, "fragments:",nrow(peaks[[pk]]$data)),...)
## }


## ##' Plot the aligned mass spec from the gatherInfo() in an averlapped way
## ##'
## ##' Plot the aligned mass spec from the gatherInfo() in an averlapped way
## ##' @title plotFragments
## ##' @name plotFragments
## ##' @param peaks object fron gatherInfo()
## ##' @param pk peak number
## ##' @param cols
## ##' @param ltys
## ##' @param TRANSFUN transorm data
## ##' @param ... further matplot() parameters
## ##' @return plots the EI mass spec
## ##' @author MR and RR
## .plotFragments <- function(peaks, pk = 1, cols = NULL, ltys = NULL, TRANSFUN = log2, ...)
## {
##     if(is.null(cols)){

##         if(ncol(peaks[[pk]]$data)/4 <= 1)
##         {
##                         cols <- c("black","green","blue","red")
##         }
##         else
##         {
##             cols <- rep(c("black","green","blue","red"), each = ncol(peaks[[pk]]$data)/4)
##         }
##     }
##     if(is.null(ltys))
##     {
##         ltys <- rep(1, length(peaks))
##     }
##     pr <- paste(round(range(peaks[[pk]]$rt, na.rm = TRUE), 2), collapse = "-")
##     matplot(peaks[[pk]]$mz, TRANSFUN(peaks[[pk]]$data), type = "h", col = cols, lty = ltys,
##             xlab = "m/z", ylab = "Intensity",
##             main = paste("index:", pk, "rt:", pr, "fragments:", nrow(peaks[[pk]]$data)), ...)
## }
