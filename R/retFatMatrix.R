#' retFatMatrix
#' 
#' Build a fat data matrix
#' 
#' This function allows to extract the data from an object created using
#' \code{gatherInfo} and build a data matrix using the area of the deconvoluted
#' and aligned peaks. The row are the samples while the column represent the
#' different peaks.
#' 
#' @param object peakDataset object
#' @param data a gatherInfo() object
#' @param minFilter the minimum number for a feature to be returned in the data
#' matrix. Default is 2/3 of the samples
#' @return A fat data matrix containing the area of the deconvoluted and
#' aligned peaks. The row are the samples while the column represent the
#' different peaks
#' @author Riccardo Romoli \email{riccardo.romoli@@unifi.it}
#' @seealso \code{\link{gatherInfo}}
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
#' data <- addXCMSPeaks(files[1:2], data, settings = mfp)
#' data
#' ma <- multipleAlignment(pd = data, group = c(1,1),
#'                         filterMin = 1, metric = 2, type = 2)
#' outList <- gatherInfo(data, ma)
#' mtxD <- retFatMatrix(object = data, data = outList, minFilter = 1)
#' 
#' @export retFatMatrix
retFatMatrix <- function (object, data, 
    minFilter = round(length(object@files) / 3 * 2)) {
        a <- lapply(seq(along = data), function(x) {
        apply(data[[x]]$data, 2, sum)
        }
        )
    ## i nomi delle colonne equivalgono al numero del file;
    ## il numero della riga equivale alla tasca della lista di gatherInfo()
    abumtx <- do.call(rbind, a)
    abumtx <- apply(abumtx, 1, "[")
    files_to_merge <- rownames(abumtx)
    if (length(grep(pattern = "^[1-9].",  files_to_merge)) == 0) {
        files.idx <- as.numeric(sub(pattern = "^.", replacement = "", 
            files_to_merge))
    }
    if (length(grep(pattern = "^[1-9].",  files_to_merge)) > 0) {
        files.idx <- as.numeric(sub(pattern = "^[1-9].", replacement = "", 
            files_to_merge))
    }
    sample <- object@files[files.idx]
    colnames(abumtx) <- sapply(1:ncol(abumtx), function(x) {
        paste0("Feat", x)
        }
        )
    mf <- minFilter
    keep <- c()
    for (g in 1:ncol(abumtx)) {
        keep[g] <- sum(!is.na(abumtx[, g])) >= mf
        }

    abumtx[is.na(abumtx)] <- c(0)
    df <- cbind.data.frame(sample, abumtx[, keep])
    return(df)
}
