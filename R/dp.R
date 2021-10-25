#' Dynamic programming algorithm, given a similarity matrix
#' 
#' This function calls C code for a bare-bones dynamic programming algorithm,
#' finding the best cost path through a similarity matrix.
#' 
#' 
#' This is a pretty standard implementation of a bare-bones dynamic programming
#' algorithm, with a single gap parameter and allowing only simple jumps
#' through the matrix (up, right or diagonal).
#' 
#' @param M similarity matrix
#' @param gap penalty for gaps
#' @param big large value used for matrix margins
#' @param verbose logical, whether to print out information
#' @return
#' 
#' \code{list} with element \code{match} with the set of pairwise matches.
#' @author Mark Robinson
#' @seealso \code{\link{normDotProduct}}
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords manip
#' @examples
#' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # read data, peak detection results
#' pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
#' pd<-addAMDISPeaks(pd,eluFiles[1:2])
#' 
#' # similarity matrix
#' r<-normDotProduct(pd@peaksdata[[1]],pd@peaksdata[[2]])
#' 
#' # dynamic-programming-based matching of peaks
#' v<-dp(r,gap=.5)
#'
#' @useDynLib flagme
#' @export dp
dp <- function(M, gap = .5, big= 10000000000, verbose = FALSE) {

  # setup score matrix
  bigr <- c(0, big * (1:ncol(M)))
  bigc <- c(big * (1:nrow(M)))
  D <- rbind(bigr, cbind(bigc, M)); # D[1, 1] <-0
  # bigr <- c(0, gap * (1:ncol(M)))
  # bigc <- c(gap * (1:nrow(M)))
  # D <- rbind(bigr, cbind(bigc, M)); # D[1, 1] <-0
  
  # setup traceback
  phi <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  phi[1, ] <- 2; phi[, 1] <- 1; phi[1, 1] <- 3
  match <- matrix(-1, nrow = nrow(M) + ncol(M), ncol = 2) # matrix to store matches

  out <- .C("dp", D = as.double(D), M = as.double(M), gap = as.double(gap), 
    phi = as.integer(phi), nr = as.integer(nrow(M)), ncol = as.integer(ncol(M)),
    match = as.integer(match), nmatch = as.integer(0), PACKAGE = "flagme")
  #list(D=matrix(out$D,ncol=ncol(D)),phi=matrix(out$phi,ncol=ncol(D)),match=1+matrix(out$match,ncol=2)[out$nmatch:1,])
  list(match = 1 + matrix(out$match, ncol = 2)[out$nmatch:1, ])
}


#' dynRT
#' 
#' Dynamic Retention Time Based Alignment algorithm, given a similarity matrix
#' 
#' This function align two chromatograms finding the maximum similarity among
#' the mass spectra
#' 
#' @param S similarity matrix
#' @return list containing the matched peaks between the two chromatograms. The
#' number represent position of the spectra in the S matrix
#' @author riccardo.romoli@@unifi.it
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
#' ## review peak picking
#' plotChrom(data, rtrange=c(7.5, 10.5), runs=c(1:2))
#' ## similarity
#' r <- ndpRT(data@peaksdata[[1]], data@peaksdata[[2]], data@peaksrt[[1]],
#'     data@peaksrt[[2]], D = 50)
#' ## dynamic retention time based alignment algorithm
#' v <- dynRT(S = r)
#' 
#' @export dynRT
dynRT <- function(S) {
    options(warn = -1)
    ## S similarity matrix
    ## move trought S to find the highest score
    S <- round(S, digits = 3)
    id <- max(table(S)) - 1
    ## filling <- which(table(S) > 200)# weak
    filling <- which(table(S) > id)
    filling.names <- as.numeric(names(filling))
    if (filling.names > 1) {
        filling.names <- 1
    }
    trace <- c()
    ## this boolean workaround is to allow to use both the ndp
    ## function and the matching function from MR and RR 
    if (filling.names == 1) {
        for (i in 1:nrow(S)) {
            trace[i] <- which(S[i, ] == min(S[i, ])) # MR
        }
    }
    if (filling.names == 0) {
        for (i in 1:nrow(S)) {
            trace[i] <- which(S[i, ] == max(S[i, ])) # RR
        }
    }
    trace.mtx <- matrix(NA, nrow = nrow(S), ncol = 2)
    trace.mtx[, 1] <- 1:nrow(S)
    trace[which(trace == 1)] <- NA
    trace.mtx[, 2] <- trace
    ## cercare i composti matchati piu di una volta e renderli univoci
    idx <- which(table(trace.mtx[, 2]) > 1) # colonne della matrice S
    names.idx <- as.numeric(names(idx)) # i doppioni
    position <- c() # i doppi con il match piu alto
    if (filling.names == 1) {
        for (k in 1:length(names.idx)) {
            position[k] <- which(S[, names.idx[k]] == min(S[, names.idx[k]]))
        }
    }
    if (filling.names == 0) {
        for (k in 1:length(names.idx)) {
            position[k] <- which(S[, names.idx[k]] == max(S[, names.idx[k]]))
        }
    }
    tutti <- which(trace %in% names.idx)
    trace.mtx[, 2][tutti[! tutti %in% position]] <- NA
    return(list(match = trace.mtx))
    options(warn = 0)
}
