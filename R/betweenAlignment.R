#' Data Structure for "between" alignment of many GCMS samples
#'
#' This function creates a "between" alignment (i.e. comparing merged peaks)
#'
#' \code{betweenAlignment} objects gives the data structure which stores the
#' result of an alignment across several "pseudo" datasets. These pseudo
#' datasets are constructed by merging the "within" alignments.
#'
#' @aliases betweenAlignment betweenAlignment-class betweenAlignment-show show,
#' betweenAlignment-method
#' @param pD a \code{peaksDataset} object
#' @param cAList \code{list} of \code{clusterAlignment} objects, one for each
#' experimental group
#' @param pAList \code{list} of \code{progressiveAlignment} objects, one for
#' each experimental group
#' @param impList \code{list} of imputation lists
#' @param filterMin minimum number of peaks within a merged peak to be kept in
#' the analysis
#' @param gap gap parameter
#' @param D retention time penalty parameter
#' @param usePeaks logical, whether to use peaks (if \code{TRUE}) or the full
#' 2D profile alignment (if \code{FALSE})
#' @param df distance from diagonal to calculate similarity
#' @param verbose logical, whether to print information
#' @param metric numeric, different algorithm to calculate the similarity
#' matrix between two mass spectrum. \code{metric=1} call
#' \code{normDotProduct()}; \code{metric=2} call \code{ndpRT()};
#' \code{metric=3} call \code{corPrt()}
#' @param type numeric, two different type of alignment function
#' @param penality penalization applied to the matching between two mass
#' spectra if \code{(t1-t2)>D}
#' @param compress logical whether to compress the similarity matrix into a
#' sparse format.
#' @return \code{betweenAlignment} object
#' @author Mark Robinson
#' @seealso \code{\link{multipleAlignment}}
#' @references Mark D Robinson (2008).  Methods for the analysis of gas
#' chromatography - mass spectrometry data \emph{PhD dissertation} University
#' of Melbourne.
#' @keywords classes
#' @examples
#'
#' 	require(gcspikelite)
#' 	## see 'multipleAlignment'
#' @importFrom stats median
#' @export betweenAlignment
betweenAlignment <- function(pD, cAList, pAList, impList, filterMin = 1,
                             gap = 0.7, D = 10, usePeaks = TRUE, df = 30,
                             verbose = TRUE, metric = 2,  type = 2,
                             penality = 0.2, compress = FALSE){
    n <- length(pAList)
    if(length(filterMin) == 1)
    {
        filterMin <- rep(filterMin, n)
    }
    pkd <- vector("list", n)
    pkr <- vector("list", n)
    filtind <- vector("list", n)  # list of filtered indices
    newind <- vector("list", n)
    for(g in 1:n){
        pa <- pAList[[g]]
        nm <- length(pa@merges)
        runs <- pa@merges[[nm]]$runs
        ind <- pa@merges[[nm]]$ind
        keep <- rowSums(!is.na(ind)) >= filterMin[g] ## RR
        if (verbose)
        {
            cat("[betweenAlignment]", names(pAList)[g], ": Removing", nrow(ind) - sum(keep),"peaks.\n")
            ind <- ind[keep,]
        }
        newind[[g]] <- lapply(impList[[g]], FUN = function(u,k) u[k,], k = keep)
        filtind[[g]] <- ind
        rt <- matrix(NA, nrow(ind), ncol(ind))
        mz <- pD@mz
        peaksdata <- matrix(0, nrow = length(pD@mz), ncol = nrow(ind))
        for(j in 1:ncol(ind)){
            if(usePeaks)
            {
                cur.ds <- pD@peaksdata[[runs[j]]]
                cur.rt <- pD@peaksrt
            }
            else
            {
                cur.ds <- pD@rawdata[[runs[j]]]
                cur.rt <- pD@rawrt
            }
            for(i in 1:nrow(ind)){
                if(!is.na(ind[i, j]))
                {
                    peaksdata[,i] <- peaksdata[,i] + cur.ds[,ind[i,j]]
                    rt[i,j] <- cur.rt[[runs[j]]][ind[i,j]]
                }
            }
        }
      pkd[[g]] <- peaksdata
	if(!usePeaks)
  {
      pkd[[g]] <- pkd[[g]] / outer(rep(1, nrow(peaksdata)), rowSums(!is.na(ind)))
  } # average over the number of samples
        pkr[[g]] <- apply(rt, 1, median, na.rm = TRUE)
    }
    # wRA <- new("peaksDataset", mz = mz, rawdata = NULL, files = names(cAList), 
    #             rawrt = NULL, peaksdata = pkd, peaksrt = pkr, filtind = filtind)
    wRA <- new("peaksDataset",
               mz = mz, rawdata = list(NULL), files = names(cAList),
               rawrt = list(NULL), peaksdata = pkd, peaksrt = pkr)
    ## return(wRA)
    ## filtind <-
    ## RR
    cA <- clusterAlignment(wRA, runs = 1:n, gap = gap, D = D, df = df,
                           metric = metric, type = type,
                           compress = compress, penality = penality) ## bug here with filterMin > 1
    pA <- progressiveAlignment(wRA, cA, gap = gap, D = D, df = df,
                               compress = compress, type = type)
    ##
    ## cA <- clusterAlignment(wRA, 1:n, gap = gap, D = D, df = df, verbose = verbose)
    ## pA <- progressiveAlignment(wRA, cA, gap = gap, D = D, df = df, verbose = verbose)

    nm <- length(pA@merges)
    ind <- pA@merges[[nm]]$ind
    groups <- pA@merges[[nm]]$runs
    full.runs <- NULL
    full.groups <- NULL
    full.ind <- matrix(NA, nrow(ind), length(pD@rawdata))
    full.newind <- vector("list", 3)
    for(i in 1:length(full.newind)){
        full.newind[[i]] <- matrix(NA, nrow(ind), length(pD@rawdata))
    }
    names(full.newind) <- c("apex", "start", "end")

    col <- 1
    for(i in 1:length(groups)){
        g <- groups[i]
        n <- length(pAList[[g]]@merges)
        full.runs <- c(full.runs, pAList[[g]]@merges[[n]]$runs)
        nr <- length(pAList[[g]]@merges[[n]]$runs)
        full.groups <- c(full.groups, rep(wRA@files[g],nr))
        this.ind <- filtind[[g]]
        rw <- which(!is.na(ind[,i]))
        cl <- col:(col + ncol(this.ind) - 1)
        full.ind[rw, cl] <- this.ind
        if(!is.null(newind[[g]]$apex))
        {
            for(j in 1:length(full.newind)){
                full.newind[[j]] <- newind[[g]][[j]]
            }
        }
        col <- col + ncol(this.ind)
    }

    new("betweenAlignment",
        mergedPeaksDataset = wRA,
        ind = full.ind,
        imputeind = full.newind,
        runs = full.runs,
        groups = full.groups,
        cA = cA,
        pA = pA,
        filtind = filtind,
        newind = newind
        )
}


#' Show method for "between" alignment of many GCMS samples
#'
#' This function show the results of a "between" alignment
#'
#' \code{betweenAlignment} objects gives the data structure which stores the
#' result of an alignment across several "pseudo" datasets. These pseudo
#' datasets are constructed by merging the "within" alignments.
#' @title betweenAlignment-show
#' @param object 
#' @return \code{betweenAlignment} object
#' @author Mark Robinson
#' @export
#' @noRd
setMethod("show","betweenAlignment",
          function(object) {
              cat("An object of class \"", class(object), "\"\n", sep = "")
              cat(length(object@mergedPeaksDataset@peaksrt), "groups:",
                  sapply(object@mergedPeaksDataset@peaksrt, length),
                  "merged peaks\n"
                  )
          }
          )
