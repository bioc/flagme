#' Add xcms/CAMERA peak detection results
#'
#' Reads the raw data using xcms, group each extracted ion according to their
#' retention time using CAMERA and attaches them to an already created
#' \code{peaksDataset} object
#'
#' Repeated calls to xcmsSet and annotate to perform peak-picking and
#' deconvolution. The peak detection results are added to the original
#' \code{peaksDataset} object. Two peak detection alorithms are available:
#' continuous wavelet transform (peakPicking=c('cwt')) and the matched filter
#' approach (peakPicking=c('mF')) described by Smith et al (2006). For further
#' information consult the xcms package manual.
#' @title addXCMSPeaks
#' @param files list of chromatogram files
#' @param object a \code{peakDataset} object
#' @param settings list. It conteins the settings for the peak-picking
#' @param rtrange vector; retention time range
#' @param mzrange vector, mz range
#' @param perfwhm etermines the maximal retentiontime difference of features in
#'     one pseudospectrum.
#' @param minintens minimum ion intensity to be included into a pseudospectra
#' @param minfeat minimum number of ion to be created a pseudospectra
#' @param multipleMatchedFilter logical Try to remove redundant peaks, in 
#' this case where there are any peaks within an absolute m/z value of 0.2 and 
#' within 3 s for any one sample in the xcmsSet (the largest peak is kept)
#' @param multipleMatchedFilterParam list. It conteins the settings for the
#' peak-picking. mz_abs represent the the mz range; rt_abs represent thert range
#' @param BPPARAM a parameter class specifying if and how parallel processing 
#' should be performed
#' @importFrom xcms xcmsRaw xcmsSet
#' @importFrom CAMERA annotate getpspectra
#' @importFrom stats aggregate
#' @export addXCMSPeaks
#' @return \code{peaksDataset} object
#' @author Riccardo Romoli \email{riccardo.romoli@@unifi.it}
#' @seealso \code{\link{peaksDataset}} \code{\link{findPeaks.matchedFilter}}
#' \code{\link{findPeaks.centWave}} \code{\link{xcmsRaw-class}}
#' @keywords manip
#' @examples
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
#'
#' @importFrom xcms xcmsRaw xcmsSet
#' @importFrom CAMERA annotate getpspectra
#' @importFrom stats aggregate
#' @export addXCMSPeaks
addXCMSPeaks <- function(files, object, settings, rtrange = NULL,
                           mzrange = NULL, perfwhm = 0.75, minintens = 100,
                           minfeat = 6, multipleMatchedFilter = FALSE, 
                           multipleMatchedFilterParam = list(
                               fwhm = c(5, 10, 20), mz_abs = 0.1, rt_abs = 3),
                           BPPARAM = bpparam()
                               ) {
    ## Rmpi tends to give many warnings that are not relevant to end
    ## users: this is an attempt to suppress this output
    owarn <- options("warn")
    on.exit(options(warn = owarn$warn))
    options(warn = -1)

    ## if an rtrange is given, we first find out which scans correspond
    ## to this and then use the scanRange argument of xcmsSet
    if (!is.null(rtrange)) {
        if (length(rtrange) != 2) {
            stop("Improper rtrange given!")
        }

        rtrange <- rtrange * 60  ## convert from minutes to seconds
        xr <- xcms::xcmsRaw(files[1])
        scanRange <- c(max(1, which(xr@scantime > rtrange[1])[1], na.rm = TRUE),
                       min(length(xr@scantime),
                           which(xr@scantime > rtrange[2])[1] - 1, na.rm = TRUE)
                           )
        allSettings <- c(list(files = files, scanrange = scanRange), settings)
    }
    else {
        allSettings <- c(list(files = files), settings)
    }

  ## peak-picking
  cat("Start peakPicking \n")
    row <- MSnbase::readMSData(files, centroided. = TRUE, mode = "onDisk",
        msLevel. = 1)
    if (class(settings)[1] == "MatchedFilterParam") {
        if (multipleMatchedFilter == TRUE) {
            settings@fwhm <- multipleMatchedFilterParam$fwhm[1]
            set1a <- xcms::findChromPeaks(row, param = settings )
            settings@fwhm <- multipleMatchedFilterParam$fwhm[2]
            set1b <- xcms::findChromPeaks(row, param = settings )
            settings@fwhm <- multipleMatchedFilterParam$fwhm[3]
            set1c <- xcms::findChromPeaks(row, param = settings )
            set1 <- set1c
            xcms::chromPeaks(set1) <- rbind(
                xcms::chromPeaks(set1a), xcms::chromPeaks(set1b), 
                xcms::chromPeaks(set1c)
            )
            xcms::chromPeaks(set1) <- xcms::chromPeaks(set1)[
                order(xcms::chromPeaks(set1)[, "sample"], decreasing = FALSE), ]
                s <- de_Duper(set1,
                             mz_abs = multipleMatchedFilterParam$mz_abs,
                             rt_abs = multipleMatchedFilterParam$rt_abs)
        }
    }
    xset <- xcms::findChromPeaks(row, param = settings, BPPARAM = bpparam())
  cat("peakPicking Done \n")
    pdp <- xcms::PeakDensityParam(sampleGroups = seq(along = files), bw = 20,
        minFraction = 0.5)
    xset <- xcms::groupChromPeaks(xset, pdp)
    xset <- xcms::fillChromPeaks(xset, param = xcms::ChromPeakAreaParam())
    xset <- as(xset, "xcmsSet")
    if (!is.null(mzrange)) {
        idx <-  (xset@peaks[, "mz"] > mzrange[1]) &
                    (xset@peaks[, "mz"] < mzrange[2]
                )
        xset@peaks <- xset@peaks[idx, ]
    }
    ## deconvolution; list of all the xset
    ap <- split(xset, factor(xcms::sampnames(xset),
        levels = xcms::sampnames(xset)))
    apd <-  lapply(ap, function(x) {
        xa <- CAMERA::xsAnnotate(x,  sample = 1)
        xa <- CAMERA::groupFWHM(xa, perfwhm = perfwhm)
    }
                    )
    ## filter pseudo spectra
    intensity <- "maxo"
    res <- lapply(apd, function(x) {
        allpks <- x@groupInfo
        minI <- minintens# * max(allpks[, intensity])
        tooSmall <- which(allpks[, intensity] < minI)
        pspectra <- lapply(x@pspectra, function(x) {
            x[!x %in% tooSmall]})
        npeaks <- sapply(pspectra, length)
        pspectra <- pspectra[npeaks >= minfeat]
        ## list of unique pspec, double masses removed
        listpspec <- lapply(pspectra, function(x) {
            aa <- cbind(mz = round(allpks[x, "mz"], digits = 0),
                        allpks[x, c(intensity, "rt", "rtmin", "rtmax")]
                        )
            double <- duplicated(aa[, 1])
            bb <- cbind(aggregate(aa[, 2], list(aa[, 1]),
                FUN = sum), aa[!double, 3:5])
            setNames(bb, c(colnames(aa)))
        }
                            )
        ## get mzrange from data
        mz.max <- max(sapply(listpspec, function(x) {
            max(x[, "mz"])}))
        mz.min <- min(sapply(listpspec, function(x) {
            min(x[, "mz"])}))
        mz.range <- data.frame(mz = c(mz.min:mz.max))
        ## merge pspec with mzrange
        listpspec.merged <- lapply(listpspec, function(x) {
            merge(x, mz.range, by = "mz", all = TRUE)
            }
            )
            }
            )
    ## merge again with a common mz range among all sampless
    max.mz <- max(sapply(res, function(x) {
        max(sapply(x, "[[", "mz"))}))
    min.mz <- min(sapply(res, function(x) {
        min(sapply(x, "[[", "mz"))}))
    mz.range.all <- data.frame(mz = c(min.mz:max.mz))
    res.mz.mrg <- lapply(res, function(x) {
        lapply(x, function(y) {
            merge(y, mz.range.all, by = "mz", all = TRUE)
            }
            )
            }
            )
    ## prepare the S4 slots
    spec.ind <- lapply(res.mz.mrg, function(x) {
        1:length(x)
        }
        )
    apex.rt <- lapply(res.mz.mrg, function(x) {
        rt <- lapply(x, "[[",  "rt")
        round(sapply(rt, function(x) {
            mean(x, na.rm = TRUE) / 60
            }
            ), digits = 3)
            }
            )
    start.rt <- lapply(res.mz.mrg, function(x) {
        rt <- lapply(x, "[[",  "rtmin")
        round(sapply(rt, function(x) {
            mean(x, na.rm = TRUE) / 60
            }
            ), digits = 3)
            }
            )
    stop.rt <- lapply(res.mz.mrg, function(x) {
        rt <- lapply(x, "[[",  "rtmax")
        round(sapply(rt, function(x) {
            mean(x, na.rm = TRUE) / 60
            }
            ), digits = 3)
            }
            )
    data <- lapply(res.mz.mrg, function(x) {
        a <- lapply(x, "[[",  intensity)
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
##' @param mz_abs mz range
##' @param rt_abs rt range
##' @return an object of xcms class
##' @author r
de_Duper <- function(object, mz_abs = 0.1, rt_abs = 2) {
    mzdiff <- 0
    peaks_mat <- xcms::chromPeaks(object)
    mz_min <- peaks_mat[, "mz"] - mz_abs
    mz_max <- peaks_mat[, "mz"] + mz_abs
    rt_min <- peaks_mat[, "rt"] - rt_abs
    rt_max <- peaks_mat[, "rt"] + rt_abs

    peaks_mat_out <- NULL

    samples <- unique(peaks_mat[, "sample"])

    cat("\n", "Duplicate peak removal; % complete: ")
    percplus <- -1

    for (i in 1:length(samples)) {
        perc <- round(i / length(samples) * 100)
        if (perc %% 10 == 0 && perc != percplus) {
            cat(perc, " ")
        }
        percplus <- perc

        peaks_mat_i <- peaks_mat[which(peaks_mat[, "sample"] == samples[i]), ,
            drop = FALSE]
        mz_min_i <- mz_min[which(peaks_mat[, "sample"] == samples[i])]
        mz_max_i <- mz_max[which(peaks_mat[, "sample"] == samples[i])]
        rt_min_i <- rt_min[which(peaks_mat[, "sample"] == samples[i])]
        rt_max_i <- rt_max[which(peaks_mat[, "sample"] == samples[i])]
        uorder_i <- order(peaks_mat_i[, "into"], decreasing = TRUE)
        uindex_i <- xcms::rectUnique(cbind(mzmin = mz_min_i, mzmax = mz_max_i,
                                            rtmin = rt_min_i, rtmax = rt_max_i),
                                      uorder_i, mzdiff)
        peaks_mat_i <- peaks_mat_i[uindex_i, , drop = FALSE]
        peaks_mat_out <- rbind(peaks_mat_out, peaks_mat_i)
    }
    cat("\n")
    xcms::chromPeaks(object) <- peaks_mat_out
    return(object)
}


#### ---------------------------------------------------------------------------
#### DEVEL
#### ---------------------------------------------------------------------------
#' Add xcms/CAMERA peak detection results
#'
#' Reads the raw data using xcms, group each extracted ion according to their
#' retention time using CAMERA and attaches them to an already created
#' \code{peaksDataset} object
#'
#' Repeated calls to xcmsSet and annotate to perform peak-picking and
#' deconvolution. The peak detection results are added to the original
#' \code{peaksDataset} object. Two peak detection alorithms are available:
#' continuous wavelet transform (peakPicking=c('cwt')) and the matched filter
#' approach (peakPicking=c('mF')) described by Smith et al (2006). For further
#' information consult the xcms package manual.
#' @title addXCMSPeaks
#' @param files list of chromatogram files
#' @param object a \code{peakDataset} object
#' @param settings list. It conteins the settings for the peak-picking
#' @param rtrange vector; retention time range
#' @param mzrange vector, mz range
#' @param perfwhm etermines the maximal retentiontime difference of features in
#'     one pseudospectrum.
#' @param minintens minimum ion intensity to be included into a pseudospectra
#' @param minfeat minimum number of ion to be created a pseudospectra
#' @importFrom xcms xcmsRaw xcmsSet
#' @importFrom CAMERA annotate getpspectra
#' @importFrom stats aggregate
#' @export addXCMSPeaks
#' @return \code{peaksDataset} object
#' @author Riccardo Romoli \email{riccardo.romoli@@unifi.it}
#' @seealso \code{\link{peaksDataset}} \code{\link{findPeaks.matchedFilter}}
#' \code{\link{findPeaks.centWave}} \code{\link{xcmsRaw-class}}
#' @keywords manip
#' @examples
#' files <- list.files(path = paste(find.package("gcspikelite"), "data",
#'                                  sep = "/"),"CDF", full = TRUE)
#' data <- peaksDataset(files[1:2], mz = seq(50, 500))
#' ## create settings object
#' mf <- list(method =  "matchedFilter", step =  0.5, steps =  2, mzdiff =  0.5,
#'            fwhm = 8, snthresh =  2, max = 500,
#'            BPPARAM = SnowParam(2, log = TRUE, stop.on.error = FALSE))
#' cwt <- list(method = "centWave", ppm = 800, snthresh = 3,
#'             mzCenterFun = "apex", peakwidth = c(3, 20),
#'             prefilter = c(1, 500),fitgauss = FALSE,
#'             firstBaselineCheck = FALSE,
#'             BPPARAM = SnowParam(2, log = TRUE, stop.on.error = FALSE))
#'
#' d1 <- addXCMSPeaks2(files[1:2], data, settings = mf, minintens = 100)
#'
addXCMSPeaks2 <- function (files, object, settings, rtrange = NULL,
                           mzrange = NULL, perfwhm = 0.75, minintens = 100,
                           minfeat = 6)
{
    ## Rmpi tends to give many warnings that are not relevant to end
    ## users: this is an attempt to suppress this output
    owarn <- options("warn")
    on.exit(options(warn = owarn$warn))
    options(warn = -1)

    ## if an rtrange is given, we first find out which scans correspond
    ## to this and then use the scanRange argument of xcmsSet
    if(!is.null(rtrange))
    {
        if(length(rtrange) != 2)
        {
            stop("Improper rtrange given!")
        }

        rtrange <- rtrange * 60  ## convert from minutes to seconds
        xr <- xcms:::xcmsRaw(files[1])
        ## range(profMz(xr)) # get mz range
        scanRange <- c(max(1,
                           which(xr@scantime > rtrange[1])[1],
                           na.rm = TRUE),
                       min(length(xr@scantime),
                           which(xr@scantime > rtrange[2])[1] - 1,
                           na.rm = TRUE))
        allSettings <- c(list(files = files, scanrange = scanRange),
                         settings)
    }
    else
    {
        allSettings <- c(list(files = files), settings)
    }

    ## peak-picking
    xset <- do.call(xcmsSet, allSettings)
    xset <- xcms:::group(xset, bw = 20, minfrac = 0.5, minsamp = 1, mzwid = 1, max = 50, sleep = 0)
    ## xset <- xcms:::fillPeaks(xset) # error
    if(!is.null(mzrange))
    {
        idx <-  (xset@peaks[,"mz"] > mzrange[1]) & (xset@peaks[,"mz"] < mzrange[2])
        xset@peaks <- xset@peaks[idx,]
    }
    ## deconvolution; list of all the xset
    ap <- split(xset, factor(sampnames(xset), levels = sampnames(xset)))
    apd <-  lapply(ap, function(x)
    {
        xa <- CAMERA:::xsAnnotate(x,  sample = 1)
        xa <- CAMERA:::groupFWHM(xa, perfwhm = perfwhm)
    })
    ## filter pseudo spectra
    intensity <- "maxo"
    res <- lapply(apd, function(x)
    {
        allpks <- x@groupInfo
        ## intensity <- "maxo"
        minI <- minintens# * max(allpks[, intensity])
        tooSmall <- which(allpks[, intensity] < minI)
        pspectra <- lapply(x@pspectra, function(x) {x[!x %in% tooSmall]})
        npeaks <- sapply(pspectra, length)
        pspectra <- pspectra[npeaks >= minfeat]
        ## list of unique pspec, double masses removed
        listpspec <- lapply(pspectra, function(x){
            aa <- cbind(mz = round(allpks[x,"mz"], digits = 0),
                        allpks[x, c(intensity, "rt", "rtmin", "rtmax")]
                        )
            double <- duplicated(aa[,1])
            bb <- cbind(aggregate(aa[, 2], list(aa[,1]), FUN = sum),aa[!double,3:5])
            setNames(bb, c(colnames(aa)))
        }
        )
        ## get mzrange from data
        mz.max <- max(sapply(listpspec, function(x){max(x[,"mz"])}))
        mz.min <- min(sapply(listpspec, function(x){min(x[,"mz"])}))
        mz.range <- data.frame(mz = c(mz.min:mz.max))
        ## merge pspec with mzrange
        listpspec.merged <- lapply(listpspec, function(x){
            merge(x, mz.range, by = "mz", all = TRUE)
        })
    })
    ## merge again with a common mz range among all sampless
    max.mz <- max(sapply(res, function(x){max(sapply(x, "[[", "mz"))}))
    min.mz <- min(sapply(res, function(x){min(sapply(x, "[[", "mz"))}))
    mz.range.all <- data.frame(mz = c(min.mz:max.mz))
    res.mz.mrg <- lapply(res, function(x){
        lapply(x, function(y){
            merge(y, mz.range.all, by = "mz", all = TRUE)
        }
        )
    })
    ## prepare the S4 slots
    spec.ind <- lapply(res.mz.mrg, function(x){1:length(x)})
    apex.rt <- lapply(res.mz.mrg, function(x){
        rt <- lapply(x, "[[",  "rt")
        round(sapply(rt, function(x){mean(x, na.rm = TRUE)/60}), digits = 3)

    })
    start.rt <- lapply(res.mz.mrg, function(x){
        rt <- lapply(x, "[[",  "rtmin")
        round(sapply(rt, function(x){mean(x, na.rm = TRUE)/60}), digits = 3)
    })
    stop.rt <- lapply(res.mz.mrg, function(x){
        rt <- lapply(x, "[[",  "rtmax")
        round(sapply(rt, function(x){mean(x, na.rm = TRUE)/60}), digits = 3)
    })
    data <- lapply(res.mz.mrg, function(x){
        a <- lapply(x, "[[",  intensity)
        aa <- do.call(cbind, a)
        colnames(aa) <- c(1:ncol(aa))
        aa[is.na(aa)] <- c(0)
        return(aa)
    })

    for(i in 1:length(files))
    {
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
deDuper <- function(object, mz.abs = 0.1, rt.abs = 2)
{
##    require("xcms")

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

    for(i in 1:length(samples))
    {
        perc <- round(i / length(samples) * 100)
        if(perc %% 10 == 0 && perc != percplus)
        {
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


#' Add xcms/CAMERA peak detection results
#'
#' Reads the raw data using xcms, group each extracted ion according to their
#' retention time using CAMERA and attaches them to an already created
#' \code{peaksDataset} object
#'
#' Repeated calls to xcmsSet and annotate to perform peak-picking and
#' deconvolution. The peak detection results are added to the original
#' \code{peaksDataset} object. Two peak detection alorithms are available:
#' continuous wavelet transform (peakPicking=c('cwt')) and the matched filter
#' approach (peakPicking=c('mF')) described by Smith et al (2006). For further
#' information consult the xcms package manual.
#'
#' @param files character vector of same length as \code{object@rawdata} (user
#' ensures the order matches)
#' @param object a \code{peaksDataset} object.
#' @param peakPicking Methods to use for peak detection. See details.
#' @param perfwhm percentage of full width half maximum. See
#' CAMERA::groupFWHM() for more details
#' @param quick logical. See CAMERA::annotate() for more details
#' @param ... arguments passed on to \code{xcmsSet} and \code{annotate}
#' @return \code{peaksDataset} object
#' @author Riccardo Romoli \email{riccardo.romoli@@unifi.it}
#' @seealso \code{\link{peaksDataset}} \code{\link{findPeaks.matchedFilter}}
#' \code{\link{findPeaks.centWave}} \code{\link{xcmsRaw-class}}
#' @keywords manip
#' @examples
#'
#' # need access to CDF (raw data)
#' require(gcspikelite)
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
#'
#' # full paths to file names
#' cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
#'
#' # create a 'peaksDataset' object and add XCMS peaks to it
#' pd <- peaksDataset(cdfFiles[1], mz=seq(50,550), rtrange=c(7.5,8.5))
#' pd <- addXCMSPeaks(cdfFiles[1], pd, peakPicking=c('mF'),
#'                    snthresh=3, fwhm=4, step=1, steps=2, mzdiff=0.5)
#'
#' @importFrom xcms xcmsRaw xcmsSet
#' @importFrom CAMERA annotate getpspectra
#' @importFrom stats aggregate
#' @export addXCMSPeaks
addXCMSPeaks <- function (files, object, peakPicking = c("cwt", "mF"),
                          multipleMatchedFilter = FALSE, perfwhm = 0.75,
                          quick = TRUE,
                          multipleMatchedFilterParam = list(fwhm = c(5, 10, 15),
                                                         mz.abs = 0.2, rt.abs = 2),
                          ...)
{
    options(warn = -1)
    cdfFiles <- as.character(files)
    if (length(cdfFiles) != length(object@rawdata))
        stop("Number of files must be the same as the number of runs (and must match).")
    xs <- lapply(cdfFiles, function(x, y) {
        f <- which(cdfFiles %in% x)
        xr <- xcmsRaw(x)
        rtrange <- c(min(object@rawrt[[f]]), max(object@rawrt[[f]])) * 60
        scanRange <- c(max(1, which(xr@scantime > rtrange[1])[1], na.rm = TRUE),
                       min(length(xr@scantime), which(xr@scantime > rtrange[2])[1] - 1, na.rm = TRUE))
        if(peakPicking == "cwt")
        {
            s <- xcmsSet(x, method = "centWave", #prefilter = c(5, 100),
                         scanrange = scanRange, integrate = 1, mzdiff = -0.001,
                         #fitgauss = TRUE,
                         ...)
        }
        if(peakPicking == "mF")
        {
            if(multipleMatchedFilter == FALSE)
            {
                s <- xcmsSet(x, method = "matchedFilter", scanrange = scanRange, max = 500, ...) # original
            }
            else
            {## make 3 xcmsSet objects using 3 FWHM values keeping all else the same
                set1a <- xcmsSet(x, method = "matchedFilter",
                                 fwhm = multipleMatchedFilterParam$fwhm[1],
                                 scanrange = scanRange, max = 500, ...)
                set1b <- xcmsSet(x, method = "matchedFilter",
                                 fwhm = multipleMatchedFilterParam$fwhm[2],
                                 scanrange = scanRange, max = 500, ...)
                set1c <- xcmsSet(x, method = "matchedFilter",
                                 fwhm = multipleMatchedFilterParam$fwhm[3],
                                 scanrange = scanRange, max = 500, ...)
                ## combine into one xcmsSet by using one of the above as a template and
                ## overriding its peaklist with a combination of all three
                set1 <- set1c
                set1@peaks <- rbind(set1a@peaks, set1b@peaks, set1c@peaks)
                set1@peaks <- set1@peaks[order(set1@peaks[, "sample"], decreasing = FALSE), ]
                ## remove redundant peaks, in this case where there are any peaks within an
                ## absolute m/z value of 0.2 and within 3 s for any one sample in the xcmsSet
                ## (the largest peak is kept)
                s <- deDuper(set1,
                             mz.abs = multipleMatchedFilterParam$mz.abs,
                             rt.abs = multipleMatchedFilterParam$rt.abs)
            }

        }
        idx <- which(s@peaks[, "mz"] > min(object@mz) & s@peaks[,
            "mz"] < max(object@mz))
        s@peaks <- s@peaks[idx, ]
        if(quick == TRUE)
        {
            a <- annotate(s, perfwhm = perfwhm, quick = quick)
        }
        if(quick == FALSE)
        {
            a <- annotate(s, perfwhm = perfwhm, cor_eic_th = 0.8,
                          pval = 0.05, graphMethod = "hcs",
                          calcIso = FALSE, calcCiS = TRUE,
                          calcCaS = FALSE)
        }
        return(a)
    }, y = peakPicking)
    if (peakPicking == "cwt") {
        area <- c("intb")
    }
    if (peakPicking == "mF") {
        area <- c("intf")
    }
    data <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
        mzrange <- object@mz
        abu <- data.frame(matrix(0, nrow = length(mzrange), ncol = length(spec.idx)))
        rownames(abu) <- mzrange
        colnames(abu) <- spec.idx
        mz <- data.frame(mz = mzrange)
        abu <- sapply(spec.idx, function(z) {
            spec <- getpspectra(xs[[x]], z)[, c("mz", area)]
            spec[, "mz"] <- round(spec[, "mz"])
            if (max(table(spec[, 1])) > 1) {
                spec.noDouble <- cbind(aggregate(spec[, 2], list(spec[,1]),
                                                 FUN = sum))
                colnames(spec.noDouble) <- c("mz", area)
                spec <- spec.noDouble
            }
            else {
                spec
            }
            abu$z <- merge(spec, mz, by = "mz", all = TRUE)[,area]
        })
        colnames(abu) <- spec.idx
        abu[is.na(abu)] <- c(0)
        return(abu)
    })
    apex.rt <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >=
            6)]
        apex.rt <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rt")]
            rt <- round(mean(spec.rt)/60, digits = 3)
        })
        return(apex.rt)
    })
    spectra.ind <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >=
            6)]
    })
    ind.start <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >=
            6)]
        rt.start <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rtmin")]
            rt <- round(mean(spec.rt), digits = 3)
        })
        return(rt.start)
    })
    ind.stop <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >=
            6)]
        rt.stop <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rtmax")]
            rt <- round(mean(spec.rt), digits = 3)
        })
        return(rt.stop)
    })
    object@files
    object@mz
    for (i in 1:length(files)) {
        ord <- order(apex.rt[[i]])
        data[[i]] <- data[[i]][, ord]
        apex.rt[[i]] <- apex.rt[[i]][ord]
        spectra.ind[[i]] <- spectra.ind[[i]][ord]
        ind.start[[i]] <- ind.start[[i]][ord]
        ind.stop[[i]] <- ind.stop[[i]][ord]
    }
    options(warn = 0)
    nm <- lapply(files, function(u) {
        sp <- strsplit(u, split = "/")[[1]]
        sp[length(sp)]
    })
    nm <- sub(".CDF$", "", nm)
    names(data) <- names(apex.rt) <- names(spectra.ind) <- names(ind.start) <- names(ind.stop) <- nm
    new("peaksDataset", files = object@files, peaksdata = data,
        peaksrt = apex.rt, peaksind = spectra.ind, peaksind.start = ind.start,
        peaksind.end = ind.stop, rawdata = object@rawdata, rawrt = object@rawrt,
        mz = object@mz)
}
