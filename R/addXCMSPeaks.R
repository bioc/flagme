##' Add xcms/CAMERA peak detection results
##'
##' Reads the raw data using xcms, group each extracted ion according to their
##' retention time using CAMERA and attaches them to an already created
##' \code{peaksDataset} object
##'
##' Repeated calls to xcmsSet and annotate to perform peak-picking and
##' deconvolution. The peak detection results are added to the original
##' \code{peaksDataset} object. Two peak detection alorithms are available:
##' continuous wavelet transform (peakPicking=c('cwt')) and the matched filter
##' approach (peakPicking=c('mF')) described by Smith et al (2006). For further
##' information consult the xcms package manual.
##' @title addXCMSPeaks
##' @param files list of chromatogram files
##' @param object a \code{peakDataset} object
##' @param settings see \code{\link[xcms]{findPeaks-matchedFilter}}
##' \code{\link[xcms]{findPeaks-centWave}}
##' @param minintens minimum ion intensity to be included into a pseudospectra
##' @param minfeat minimum number of ion to be created a pseudospectra
##' @param BPPARAM a parameter class specifying if and how parallel processing
##' should be performed
##' @param multipleMF logical Try to remove redundant peaks, in
##' this case where there are any peaks within an absolute m/z value of 0.2 and
##' within 3 s for any one sample in the xcmsSet (the largest peak is kept)
##' @param multipleMFParam list. It conteins the settings for the
##' peak-picking. mz_abs represent the the mz range; rt_abs represent thert range
##' @param mz.abs mz range
##' @param rt.abs rt range
##' @return \code{peaksDataset} object
##' @importFrom xcms xcmsRaw xcmsSet
##' @importFrom CAMERA annotate getpspectra
##' @importFrom stats aggregate
##' @export addXCMSPeaks
##' @author Riccardo Romoli \email{riccardo.romoli@unifi.it}
##' @seealso \code{\link[xcms]{peaksDataset}} \code{\link[xcms]{findPeaks-matchedFilter}}
##' \code{\link[xcms]{findPeaks-centWave}} \code{\link[xcms]{xcmsRaw-class}}
##' @keywords manip
##' @examples
##' files <- list.files(path = paste(find.package("gcspikelite"), "data",
##'                     sep = "/"),"CDF", full = TRUE)
##' data <- peaksDataset(files[1:2], mz = seq(50, 550), rtrange = c(7.5, 8.5))
##' ## create settings object
##' mfp <- xcms::MatchedFilterParam(fwhm = 10, snthresh = 5)
##' cwt <- xcms::CentWaveParam()
##' data <- addXCMSPeaks(files[1:2], data, settings = mfp, multipleMF = FALSE)
##' data
addXCMSPeaks <- function(files, object, settings = list(), minintens = 100,
                         minfeat = 6, BPPARAM = bpparam(),
                         multipleMF = FALSE, multipleMFParam = list(
                                               fwhm = c(5, 10, 15),
                                               mz.abs = 0.2, rt.abs = 2)) {
  options(warn = -1)
  cdfFiles <- as.character(files)
  if (length(cdfFiles) != length(object@rawdata))
    stop("Number of files must be the same as the number of runs
(and must match).")
  if (class(settings)[1] == "CentWaveParam") {
    all_settings <- list(files = files, method = "centWave",
                         ppm = settings@ppm, peakwidth = settings@peakwidth,
                         snthresh = settings@snthresh,
                         prefilter = settings@prefilter,
                         mzCenterFun = settings@mzCenterFun,
                         integrate = settings@integrate,
                         mzdiff = settings@mzdiff,
                         fitgauss = settings@fitgauss, noise = settings@noise,
                         firstBaselineCheck = settings@firstBaselineCheck,
                         roiScales = settings@roiScales, BPPARAM = BPPARAM
                         )
    xs <- do.call(xcms::xcmsSet, all_settings)
  } else if (class(settings)[1] == "MatchedFilterParam" && multipleMF == FALSE) {
    all_settings <- list(files = files, method = "matchedFilter",
                         fwhm = settings@fwhm, BPPARAM = BPPARAM,
                         sigma = settings@sigma, max = settings@max,
                         snthresh = settings@snthresh, steps = settings@steps,
                         mzdiff = settings@mzdiff, index = settings@index
                         )
    xs <- do.call(xcms::xcmsSet, all_settings)
  } else if (class(settings)[1] == "MatchedFilterParam" && multipleMF == TRUE) {
    all_settings <- list(files = files, method = "matchedFilter",
                         fwhm = settings@fwhm,
                         sigma = settings@sigma, max = settings@max,
                         snthresh = settings@snthresh, steps = settings@steps,
                         mzdiff = settings@mzdiff, index = settings@index
                         )
    ## multipleMF
    ## make 3 xcmsSet objects using 3 FWHM values keeping all else the same
    all_settings$fwhm <- multipleMFParam$fwhm[1]
    set1a <- do.call(xcms::xcmsSet, all_settings)
    all_settings$fwhm <- multipleMFParam$fwhm[2]
    set1b <- do.call(xcms::xcmsSet, all_settings)
    all_settings$fwhm <- multipleMFParam$fwhm[3]
    set1c <- do.call(xcms::xcmsSet, all_settings)
    ## combine into one xcmsSet by using one of the above as a template and
    ## overriding its peaklist with a combination of all three
    set1 <- set1c
    set1@peaks <- rbind(set1a@peaks, set1b@peaks, set1c@peaks)
    set1@peaks <- set1@peaks[order(set1@peaks[, "sample"], decreasing = FALSE), ]
    ## remove redundant peaks, in this case where there are any peaks within an
    ## absolute m/z value of 0.2 and within 3 s for any one sample in the xcmsSet
    ## (the largest peak is kept)
    xs <- deDuper(set1, mz.abs = multipleMFParam$mz.abs,
                  rt.abs = multipleMFParam$rt.abs)
  }
  ## xs <- do.call(xcms::xcmsSet, all_settings)
  xs <- xcms::group(xs, bw = 20, minfrac = 0.5, minsamp = 1, mzwid = 1,
                     max = 500, sleep = 0)
  xs <- xcms::fillPeaks(xs)

  ## if(!is.null(mzrange)) {
  ##   idx <-  (xs@peaks[,"mz"] > mzrange[1]) & (xs@peaks[,"mz"] < mzrange[2])
  ##   xs@peaks <- xs@peaks[idx,]
  ## }

  ## deconvolution; list of all the xset
  ap <- split(xs, factor(sampnames(xs), levels = sampnames(xs)))
  apd <-  lapply(ap, function(x) {
    xa <- CAMERA::xsAnnotate(x)
    xa <- CAMERA::groupFWHM(xa, perfwhm = 0.6)
    xa <- CAMERA::groupCorr(xa, cor_eic_th = 0.9, calcCiS = TRUE)
  }
  )

  ## filter pseudo spectra
  intensity <- "maxo"
  res <- lapply(apd, function(x) {
    allpks <- x@groupInfo
    ## intensity <- "maxo"
    minI <- minintens# * max(allpks[, intensity])
    tooSmall <- which(allpks[, intensity] < minI)
    pspectra <- lapply(x@pspectra, function(x) {x[!x %in% tooSmall]})
    npeaks <- sapply(pspectra, length)
    pspectra <- pspectra[npeaks >= minfeat]
    ## list of unique pspec, double masses removed
    listpspec <- lapply(pspectra, function(x) {
      aa <- cbind(mz = round(allpks[x,"mz"], digits = 0),
                  allpks[x, c(intensity, "rt", "rtmin", "rtmax")]
                  )
      double <- duplicated(aa[, 1])
      bb <- cbind(aggregate(aa[, 2], list(aa[,1]), FUN = sum), aa[!double, 3:5])
      setNames(bb, c(colnames(aa)))
    }
    )
    ## get mzrange from data
    mz.max <- max(sapply(listpspec, function(x) {max(x[, "mz"])}))
    mz.min <- min(sapply(listpspec, function(x) {min(x[, "mz"])}))
    mz.range <- data.frame(mz = c(mz.min:mz.max))
    ## merge pspec with mzrange
    listpspec.merged <- lapply(listpspec, function(x) {
      merge(x, mz.range, by = "mz", all = TRUE)
    }
    )
  }
  )

  ## merge again with a common mz range among all sampless
  max.mz <- max(sapply(res, function(x) {max(sapply(x, "[[", "mz"))}))
  min.mz <- min(sapply(res, function(x) {min(sapply(x, "[[", "mz"))}))
  mz.range.all <- data.frame(mz = c(min.mz:max.mz))
  res.mz.mrg <- lapply(res, function(x) {
    lapply(x, function(y) {
      merge(y, mz.range.all, by = "mz", all = TRUE)
    }
    )
  }
  )

  ## prepare the S4 slots
  spec.ind <- lapply(res.mz.mrg, function(x){1:length(x)})
  apex.rt <- lapply(res.mz.mrg, function(x){
    rt <- lapply(x, "[[", "rt")
    round(sapply(rt, function(x) {mean(x, na.rm = TRUE) / 60}), digits = 3)
  }
  )
  start.rt <- lapply(res.mz.mrg, function(x) {
    rt <- lapply(x, "[[", "rtmin")
    round(sapply(rt, function(x) {mean(x, na.rm = TRUE) / 60}), digits = 3)
  }
  )
  stop.rt <- lapply(res.mz.mrg, function(x) {
    rt <- lapply(x, "[[", "rtmax")
    round(sapply(rt, function(x) {mean(x, na.rm = TRUE) / 60}), digits = 3)
  }
  )
  data <- lapply(res.mz.mrg, function(x) {
    a <- lapply(x, "[[", intensity)
    aa <- do.call(cbind, a)
    colnames(aa) <- c(1:ncol(aa))
    aa[is.na(aa)] <- c(0)
    return(aa)
  }
  )

  for (i in 1:length(files)) {
    ord <- order(apex.rt[[i]])
    data[[i]] <- data[[i]][, ord]
    apex.rt[[i]] <- apex.rt[[i]][ord]
    spec.ind[[i]] <- spec.ind[[i]][ord]
    start.rt[[i]] <- start.rt[[i]][ord]
    stop.rt[[i]] <- stop.rt[[i]][ord]
  }

  new("peaksDataset",
      files = files,
      peaksdata = data,
      peaksrt = apex.rt,
      peaksind = spec.ind,
      peaksind.start = start.rt,
      peaksind.end = stop.rt,
      rawdata = object@rawdata,
      rawrt = object@rawrt,
      mz = c(min.mz:max.mz)
      )
}

##' Duplicate peak removal function
##'
##' Remove redundant peaks, in this case where there are any peaks within an
##' absolute m/z value of 0.2 and within 3 s for any one sample in the xcmsSet
##' (the largest peak is kept)
##' @title deDuper
##' @param object xcms object
##' @param mz.abs mz range
##' @param rt.abs rt range
##' @return an object of xcms class
##' @author r
deDuper <- function(object, mz.abs = 0.1, rt.abs = 2) {
  mzdiff = 0
  peaks.mat <- object@peaks
  mz.min <- peaks.mat[, "mz"] - mz.abs
  mz.max <- peaks.mat[, "mz"] + mz.abs
  rt.min <- peaks.mat[, "rt"] - rt.abs
  rt.max <- peaks.mat[, "rt"] + rt.abs
  peaks.mat.out <- NULL
  samples <- unique(peaks.mat[,"sample"])
  cat("\n", "Duplicate peak removal; % complete: ")
  percplus <- -1

  for (i in 1:length(samples)) {
    perc <- round(i / length(samples) * 100)
    if (perc %% 10 == 0 && perc != percplus) {
      cat(perc, " ")
    }

    percplus <- perc
    peaks.mat.i <- peaks.mat[which(peaks.mat[, "sample"] == samples[i]), , drop = FALSE]
    mz.min.i <- mz.min[which(peaks.mat[, "sample"] == samples[i])]
    mz.max.i <- mz.max[which(peaks.mat[, "sample"] == samples[i])]
    rt.min.i <- rt.min[which(peaks.mat[, "sample"] == samples[i])]
    rt.max.i <- rt.max[which(peaks.mat[, "sample"] == samples[i])]

    uorder.i <- order(peaks.mat.i[, "into"], decreasing = TRUE)
    uindex.i <- xcms:::rectUnique(cbind(mzmin = mz.min.i, mzmax = mz.max.i,
                                        rtmin = rt.min.i, rtmax = rt.max.i),
                                  uorder.i, mzdiff)
    peaks.mat.i <- peaks.mat.i[uindex.i, , drop = FALSE]
    peaks.mat.out <- rbind(peaks.mat.out, peaks.mat.i)
  }

  cat("\n")
  object@peaks <- peaks.mat.out
  return(object)
}
