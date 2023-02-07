##' Decompress method for progressiveAlignment
##'
##' Decompress method for progressiveAlignment
##' @title Compress method for progressiveAlignment
##' @param object progressiveAlignment object
##' @param verbose logical
##' @param ... dummy
##' @author MR
##' @importFrom SparseM as.matrix
##' @import methods
##' @importFrom methods setMethod new
##' @keywords internal
setMethod("decompress", "progressiveAlignment",
          function(object, verbose = TRUE, ...) {
              if (object@merges[[1]]$compressed == FALSE) {
                  if (verbose)
                      cat("[decompress.progressiveAlignment] 
                      Already decompressed.\n")
                  return(object)
              }
              for (i in 1:length(object@merges)) {
                  if (object@merges[[i]]$compressed) {
                      object@merges[[i]]$r <- 
                        1 - as.matrix(object@merges[[i]]$r)
                      object@merges[[i]]$compressed <- FALSE
                  }
              }
              new("progressiveAlignment", object)
          }
          )


##' Decompress method for progressiveAlignment
##'
##' Deompress method for progressiveAlignment
##' @title Compress method for progressiveAlignment
##' @param object dummy
##' @param verbose dummy
##' @param ... dummy
##' @author MR
##' @importFrom SparseM as.matrix.csc
##' @import methods
##' @importFrom methods setMethod new
##' @keywords internal
setMethod("compress", "progressiveAlignment",
          function(object, verbose = TRUE, ...) {
              if (object@merges[[1]]$compressed) {
                  if (verbose)
                      cat("[compress.progressiveAlignment]
                      Already compressed.\n")
                  return(object)
              }
              for (i in 1:length(object@merges)) {
                  if (!(object@merges[[i]]$compressed)) {
                      object@merges[[i]]$r <-
                        as.matrix.csc(1 - object@merges[[i]]$r)
                      object@merges[[i]]$compressed <- TRUE
                  }
              }
              new("progressiveAlignment", object)
          }
          )


##' Show method for progressiveAlignment object
##'
##' Show method for progressiveAlignment object
##' @param object progressiveAlignment object
##' @author MR
##' @export
##' @noRd
setMethod("show", "progressiveAlignment", function(object){
   cat("An object of class \"", class(object), "\"\n", sep = "")
   cat(length(object@merges), "merges\n")
}
)



#' Data Structure for progressive alignment of many GCMS samples
#'
#' Performs a progressive peak alignment (clustalw style) of multiple GCMS peak
#' lists
#'
#' The progressive peak alignment we implemented here for multiple GCMS peak
#' lists is analogous to how \code{clustalw} takes a set of pairwise sequence
#' alignments and progressively builds a multiple alignment.  More details can
#' be found in the reference below.
#'
#' @name progressiveAlignment-class
#' @aliases progressiveAlignment-class progressiveAlignment-show
#' progressiveAlignment show,progressiveAlignment-method
#' @param pD a \code{peaksDataset} object
#' @param cA a \code{clusterAlignment} object
#' @param D retention time penalty
#' @param gap gap parameter
#' @param verbose logical, whether to print information
#' @param usePeaks logical, whether to use peaks (if \code{TRUE}) or the full
#' 2D profile alignment (if \code{FALSE})
#' @param df distance from diagonal to calculate similarity
#' @param compress logical, whether to store the similarity matrices in sparse
#' form
#' @param type numeric, two different type of alignment function
#' @return \code{progressiveAlignment} object
#' @author Mark Robinson
#' @seealso \code{\link{peaksDataset}}, \code{\link{multipleAlignment}}
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
#' data <- addXCMSPeaks(files[1:2], data, settings = mfp)
#' data
#' ca <- clusterAlignment(data, gap = 0.5, D = 0.05, df = 30, metric = 1,
#'   type = 1, compress = FALSE)
#' pa <- progressiveAlignment(data, ca, gap = 0.6, D = 0.1, df = 30,
#'  type = 1, compress = FALSE)
#'
#' @export
progressiveAlignment <- function(pD, cA, D = 50, gap = 0.5, verbose = TRUE,
                                 usePeaks = TRUE, df = 30, compress = FALSE,
                                 type = 2) {
    m <- cA@merge
    merges <- vector("list", nrow(m))
    if (usePeaks) {
      pd <- pD@peaksdata
    } else {
      pd <- pD@rawdata
    }
    pD <- NULL #?
    cA <- decompress(cA, verbose = verbose)

    for (i in 1:nrow(m)) {
      if (verbose) {
            cat("[progressiveAlignment] Doing merge", m[i, ], "\n")
        }
      ## left
      if (m[i, 1] < 0) {
        left.runs <- (-m[i, 1])
        left.ind <- matrix(1:ncol(pd[[left.runs]]), ncol = 1)
        } else {
          left.runs <- merges[[m[i, 1]]]$runs
          left.ind <- merges[[m[i, 1]]]$ind
        }
      ## right
      if (m[i, 2] < 0) {
        right.runs <- abs(m[i, 2])
        # subscript out of bound
        right.ind <- matrix(1:ncol(pd[[right.runs]]), ncol = 1)
        } else {
          right.ind <- merges[[m[i, 2]]]$ind
          right.runs <- merges[[m[i, 2]]]$runs
          }
      if (verbose) {
        cat("[progressiveAlignment] left.runs:", left.runs, ", ")
        cat("right.runs:", right.runs, "\n")
      }
      if (length(right.runs) == 1 & length(left.runs) == 1) {
      al <- cA@alignments[[cA@aligned[left.runs, right.runs]]]
      v <- al@v
      sim <- al@r
      mi <- .merge.indices(nrow(left.ind), nrow(right.ind), v$match)
      new.ind <- mi
      } else {
        sim <- matrix(0, nrow = nrow(left.ind), ncol = nrow(right.ind))
        count <- matrix(0, nrow = nrow(left.ind), ncol = nrow(right.ind))
        if (verbose) {
          cat("[progressiveAlignment] (dot=50) going to", nrow(sim), ":")
          }
          for (r in 1:nrow(sim)) {
            if (verbose) {
              if(r %% 50 == 0) {
                cat(".")
                }
                }
                for (cc in max(1, r - df) : min(r + df, ncol(sim))) {
                  # cc generate subsrcipt out of bound
                  if (cc > dim(sim)[2]) {
                    cat("\n\n", "Try to increase the df parameter to get a 
                    better alignment", "\n\n")
              # cc <- dim(sim)[2]
              }
              for (lr in 1:length(left.runs)) {
                for(rr in 1:length(right.runs)) {
                  ai <- cA@aligned[left.runs[lr], right.runs[rr]]
                  # this.sim <- cA$alignments[[ai]]$r
                  lind <- left.ind[r, lr]
                  rind <- right.ind[cc, rr]
                  # Error: subscript out of bounds. Solution: increase df
                  if (!is.na(lind) & !is.na(rind)) {
                    if(left.runs[lr] < right.runs[rr]) {
                      sim[r, cc] <- sim[r, cc] + cA@alignments[[ai]]@r[lind, rind]
                      } else {
                        sim[r, cc] <- sim[r, cc] + cA@alignments[[ai]]@r[rind, lind]
                      }
                      count[r, cc] = count[r, cc] + 1
                      }
                  }
                }
              }
            }
            if (verbose) {
              cat("\n")
              } 
              if (type == 2) {
                v <- dynRT(S = sim)
                # remove non-matched peaks
                v$match <- v$match[!is.na(v$match[, 2]), ]
                }
              if (type == 1) {
                sim[count > 0] <- sim[count > 0] / count[count > 0]
                sim[sim == 0] <- 1
                v <- dp(sim, gap = gap, verbose = verbose)
                }
    mi <- .merge.indices(nrow(left.ind), nrow(right.ind), v$match)
    new.ind <- cbind(left.ind[mi[, 1], ], right.ind[mi[, 2], ])
    }
    rownames(new.ind) <- 1:nrow(new.ind)
    merges[[i]] <- list(ind = new.ind, mi = mi, runs = c(left.runs, right.runs),
      left = left.runs, right = right.runs, r = sim, compressed = FALSE)
      }
  cA <- NULL
  if (verbose) {
      print(gc())
    }
  pA <- new("progressiveAlignment", merges = merges)
  if (compress) {
      return(compress(pA, verbose = verbose))
    } else {
      return(pA)
    }
}


.merge.indices <- function(nl, nr, m) {
  lind <- cbind(1:nl, 0)
  lind[lind[, 1] %in% m[, 1], 2] <- 1
  rind <- cbind(1:nr, 0)
  rind[rind[, 1] %in% m[, 2], 2] <- 1
  mg <- matrix(NA, nrow = nrow(lind) + nrow(rind) - nrow(m), ncol = 2)
  li <- 1; ri <- 1; i <- 1
  while (i <= nrow(mg)) {
    if (li > nrow(lind) & ri <= nrow(rind)) {
      mg[i, ] <- c(NA, rind[ri, 1])
      ri <- ri + 1
      i <- i + 1
      next
      }
    if (li <= nrow(lind) & ri > nrow(rind)) {
      mg[i, ] <- c(lind[li, 1], NA)
      li <- li + 1
      i <- i + 1
      next
      }
    if (lind[li, 2] == 1) {
      if(rind[ri, 2] == 1) { # match
        mg[i, ] <- c(lind[li, 1], rind[ri, 1])
        ri <- ri + 1
		    li <- li + 1
        # cat("match", i, li, ri, "-", mg[i, ], "\n")
        } else { # right unmatched
          mg[i, ] <- c(NA, rind[ri, 1])
          ri <- ri + 1
		      # cat("right unmatched", i, li, ri, "-", mg[i, ], "\n")
          }
      } else {
        if (rind[ri, 2] == 1) { # left unmatched
          mg[i, ] <- c(lind[li, 1], NA)
		      li <- li + 1
          # cat("left unmatched", i, li, ri, "-", mg[i, ], "\n")
          } else {# both unmatched
            mg[i, ] <- c(lind[li, 1], NA)
            # cat("both unmatched - A", i, li, ri, "-", mg[i, ], "\n")
            i <- i + 1
            mg[i, ] <- c(NA, rind[ri, 1])
            li <- li + 1
            ri <- ri + 1
            # cat("both unmatched - B", i, li, ri, "-", mg[i, ], "\n")
            }
        }
        i <- i + 1
    }
    mg
}
