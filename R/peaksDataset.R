#' Data Structure for raw GCMS data and peak detection results
#'
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#'
#'
#' peaksDataset is a hold-all data structure of the raw and peak detection
#' data.
#'
#' @title peaksDataset-show
#' @param object an object of class \code{\link{peakDataset}}
#' @return \code{peaksDataset} object
#' @author Mark Robinson
#' @export
#' @noRd
setMethod("show","peaksDataset",
function(object) {
        cat("An object of class \"",class(object),"\"\n",sep="")
        cat(length(object@rawdata), "samples:",names(object@rawdata),"\n",sep=" ")
        cat(length(object@mz), "m/z bins - range: (",range(object@mz),")\n",sep=" ")
        cat("scans:",sapply(object@rawdata,ncol),"\n",sep=" ")
        cat("peaks:",sapply(object@peaksdata,ncol),"\n",sep=" ")
}) 


#' Data Structure for raw GCMS data and peak detection results
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' 
#' peaksDataset is a hold-all data structure of the raw and peak detection
#' data.
#' 
#' @aliases peaksDataset-class peaksDataset-show peaksDataset-plot peaksDataset
#' show,peaksDataset-method plot,peaksDataset-method
#' plot,peaksDataset,ANY-method
#' @param fns character vector, filenames of raw data in CDF format.
#' @param verbose logical, if \code{TRUE} then iteration progress information
#' is output.
#' @param mz vector giving bins of raw data table.
#' @param rtDivide number giving the amount to divide the retention times by.
#' @param rtrange retention time range to limit data to (must be \code{numeric}
#' vector of length 2)
#' @return
#' 
#' \code{peaksDataset} object
#' @author Mark Robinson
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords classes
#' @examples
#' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # read data
#' pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
#' show(pd)
#' 
#' @export peaksDataset
peaksDataset <- function(fns=dir(, "[Cc][Dd][Ff]"), verbose=TRUE,
                            mz=seq(50, 550), rtDivide=60, rtrange=NULL){
    rawdata <- vector("list", length(fns))
    rawrt <- vector("list", length(fns))
    for(i in 1:length(fns)){
        if(verbose) 
            cat(" Reading ", fns[i], "\n")
        a <- xcmsRaw(fns[i])
        D1 <- cbind.data.frame(mz=a@mzrange[1]:a@mzrange[2],
                               a@env$profile)
        D2 <- data.frame(mz=mz)
        if(a@mzrange[1]-range(mz)[1] < 3 & a@mzrange[2]-range(mz)[2] <
           3){
            mrg <- merge(D1, D2, by='mz', all=TRUE)
        }else{
            mrg <- merge(D1, D2, by='mz', all=FALSE)
        }# not a very clean approach. The problem is in the raw m/z;
         # the possible cause is the centroid approx. diring the
         # acquisition of the cgromatogram se all=TRUE funziona con
         # mzrange= 50-550; se all=FALSE funziona con mzrange= 100-400 
        mrg[is.na(mrg)] <- 0
        rownames(mrg) <- 1:nrow(mrg)
        colnames(mrg) <- 1:ncol(mrg)
        a@env$profile <- data.matrix(mrg)
        a@mzrange <- range(mz) ## allows mzrange param without errors
        if(is.null(mz)) 
            mz <- seq(a@mzrange[1], a@mzrange[2])
        this.mz <- seq(a@mzrange[1], a@mzrange[2])
        rawrt[[i]] <- a@scantime/rtDivide
        nc <- length(rawrt[[i]])
        if(length(rtrange) == 2){
            w <- which(rawrt[[i]] >= rtrange[1] & rawrt[[i]] <= rtrange[2])
            rawrt[[i]] <- rawrt[[i]][w]
        }else{
            w <- seq(1, nc)
        }
        rawdata[[i]] <- a@env$profile[, w]
        rm(a)
        if(length(this.mz) != length(mz)){
            d <- matrix(0, nrow=length(mz), ncol=length(w))
            d[match(this.mz, mz), ] <- rawdata[[i]]
            rawdata[[i]] <- d
        }
    }
    gc()
    nm <- lapply(fns, function(u){
        sp <- strsplit(u, split="/")[[1]]
        sp[length(sp)]
    })
    nm <- sub(".CDF", "", nm)
    names(rawdata) <- names(rawrt) <- nm
    new("peaksDataset", rawdata=rawdata,
        rawrt=rawrt, mz=mz, files=fns)
    ## new("peaksDataset",
    ##   rawdata = lapply(rawdata, function(x){ as.matrix.csc(x) }),
    ##   rawrt = rawrt,
    ##   mz = mz,
    ##   files = fns)

}




#' Plotting functions for GCMS data objects
#' 
#' Store the raw data and optionally, information regarding signal peaks for a
#' number of GCMS runs
#' 
#' Each TIC is scale to the maximum value (as specified by the
#' \code{how.near} and \code{max.near} values). The many parameters gives
#' considerable flexibility of how the TICs can be visualized.
#' 
##' @param object a \code{peaksDataset} object.
##' @param runs set of run indices to plot
##' @param mzind set of mass-to-charge indices to sum over (default, all)
##' @param mind matrix of aligned indices
##' @param plotSampleLabels logical, whether to display sample labels
##' @param calcGlobalMax logical, whether to calculate an overall maximum
##' for scaling
##' @param peakCex character expansion factor for peak labels
##' @param plotPeaks logical, whether to plot hashes for each peak
##' @param plotPeakBoundaries logical, whether to display peak boundaries
##' @param plotPeakLabels logical, whether to display peak labels
##' @param plotMergedPeakLabels logical, whether to display 'merged' peak
##' labels
##' @param mlwd line width of lines indicating the alignment
##' @param usePeaks logical, whether to plot  alignment of peaks
##' (otherwise, scans)
##' @param plotAcrossRuns logical, whether to plot across peaks when
##' unmatched peak is given
##' @param overlap logical, whether to plot TIC/XICs overlapping
##' @param rtrange  vector of length 2 giving start  and end of the X-axis
##' @param cols  vector of colours (same length as the length of runs)
##' @param thin  when \code{usePeaks=FALSE}, plot  the alignment lines
##' every \code{thin} values
##' @param max.near where to look for maximum
##' @param how.near how far away from \code{max.near} to look
##' @param scale.up  a constant factor to scale the TICs
##' @param ... further arguments passed to the \code{plot}
##' @return plot the chromatograms
##' @author Mark Robinson
##' @seealso \code{\link{peaksDataset}}
##' @references Mark D Robinson (2008).  Methods for the analysis of gas
##' chromatography - mass spectrometry data \emph{PhD dissertation} University
##' of Melbourne.
##' @keywords classes 
##' @importFrom graphics lines axis plot points text par
##' @examples
##' 
##' require(gcspikelite)
##' 
##' ## paths and files
##' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
##' cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
##' eluFiles <- dir(gcmsPath, "ELU", full=TRUE)
##' 
##' ## read data
##' pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,8.5))
##' 
##' ## image plot
##' plotChrom(pd, rtrange = c(7.5,8.5), plotPeaks = TRUE,
##'           plotPeakLabels = TRUE)
##' @export
setMethod("plotChrom","peaksDataset",
          function(object, runs = 1:length(object@rawdata),
                   mzind = 1:nrow(object@rawdata[[1]]),
                   mind = NULL, plotSampleLabels = TRUE,
                   calcGlobalMax = FALSE, peakCex = 0.8,
                   plotPeaks = TRUE, plotPeakBoundaries = FALSE,
                   plotPeakLabels = FALSE,
                   plotMergedPeakLabels = TRUE, mlwd = 3,
                   usePeaks = TRUE, plotAcrossRuns = FALSE,
                   overlap = F, rtrange = NULL, cols = NULL, thin = 1,
                   max.near = median(object@rawrt[[1]]), how.near = 50,
                   scale.up = 1, ...){
    ## only specific one of ma or runs
  if(is.null(rtrange))
    rtrange <- c(min(object@rawrt[[1]]), max(object@rawrt[[1]]))
  if(is.factor(cols))
    cols <- c("black", "green", "blue", "red", "grey", "orange")[cols]
  if(is.null(cols))
      cols <- rep(c("black", "green", "blue", "red"),
                  each = length(object@rawdata)/4)
  if(length(cols) != length(runs))
    cols <- rep("black", length(runs))
  if(usePeaks){
      datalist <- object@peaksdata
      rtlist <- object@peaksrt
  }else{
      datalist <- object@rawdata
      rtlist <- object@rawrt
  }
  if(calcGlobalMax){
      mx <- 0
      for(i in 1:length(runs)) {
          ind <- runs[i]
          mxind <- which(abs(object@rawrt[[ind]]-max.near)<how.near)
          b <- colSums(matrix(object@rawdata[[ind]][mzind,],
                              nrow = length(mzind)))
          if(mx < max(b[mxind]))
              mx <- max(b[mxind])
      }
  }
  if(length(object@peaksdata) != length(object@rawdata))
    plotPeaks <- FALSE
  for(i in 1:length(runs)){
      ind <- runs[i]
      b <- colSums(matrix(object@rawdata[[ind]][mzind,],
                          nrow = length(mzind)))
      if(!calcGlobalMax){
          mxind <- which(abs(object@rawrt[[ind]]-max.near) < how.near)
          mx <- max(b[mxind])
      }
    b <- b / mx*scale.up
	if(!overlap){
      if(i == 1){
          plot(object@rawrt[[ind]], b, type = "l", xlab = "Retention Time",
               ylab = "", xlim = rtrange,
               ylim = c(-.5,length(runs)), col = cols[i], xaxs = "i",
               axes = FALSE, ...)
            axis(1)
	  }else{
	    lines(object@rawrt[[ind]], b+i-1, col = cols[i], ...)
	  }
	  if(plotSampleLabels)
        text(rtrange[1]+.1*diff(rtrange), i-.5, names(object@rawdata)[ind],
             col = "black")
    } else {
      if(i == 1) {
          plot(object@rawrt[[ind]], b, type = "l", xlab = "Retention Time",
               ylab = "", xlim = rtrange,
               ylim = c(0,1), col = cols[i], xaxs = "i", axes = FALSE)
            axis(1)
	  } else {
	    lines(object@rawrt[[ind]], b, col = cols[i])
	  }
	}
    if(plotPeaks & !overlap){
   	  rt <- rtlist[[ind]]; d <- min(diff(rt))*.4
	  for(j in 1:length(rt)) {
	    if(rt[j] > rtrange[1] & rt[j] < rtrange[2]) {
		  if (plotPeakBoundaries) {
			st <- object@rawrt[[ind]][ object@peaksind.start[[ind]][j] ]
			en <- object@rawrt[[ind]][ object@peaksind.end[[ind]][j] ]
      ## cat(ind,st,en,rt[j],"\n")
			if(j%%2==0)
      {
          lines(c(st,st,rt[j],rt[j],rt[j],en,en),
                c(i-1+.025,i-1,i-1,i-1+.3,i-1,i-1,i-1+.025), col = "blue")
      }
			else
          {
              lines(c(st,st,rt[j],rt[j],rt[j],en,en),
                c(i-1+.025,i-1,i-1,i-1+.3,i-1,i-1,i-1+.025)-.025, col = "red")
                ## lines(c(rt[j],rt[j]),c(i-1,i-1+.3))
          }
		  } else {
	        lines(c(rt[j]-d,rt[j]+d), c(i-1,i-1))
	        lines(c(rt[j],rt[j]), c(i-1,i-1+.3))
		  }
	    }
	  }
	  if(plotPeakLabels)
	    text(rt, i-1+.15, labels = 1:length(rt), adj = 0, cex = peakCex)
    }
  }
  if(!is.null(mind) & !overlap) {
	  rws <- as.integer(seq(1, nrow(mind), by = thin))
	  for(i in rws) {
	    if( !is.na(mind[i,1]) & plotMergedPeakLabels )
          text(rtlist[[runs[1]]][mind[i,1]], 0, srt = 90, i, col = "blue",
               pos = 1, cex = peakCex*.8, adj = 1)
	    for(j in 1:(ncol(mind)-1)) {
	      ind1 <- mind[i,j]
	      ind2 <- mind[i,j+1]
		  if(!is.na(ind1) & !is.na(ind2)) {
		    rt1 <- rtlist[[runs[j]]][ind1]
        rt2 <- rtlist[[runs[j+1]]][ind2]
			## cat(ind1,ind2,runs[j],runs[j+1],rt1,rtrange[1],rt2,rtrange[2],"\n")
        if((rt1 > rtrange[1] | rt2 > rtrange[1]) &
           (rt1 < rtrange[2] | rt2 < rtrange[2])){
			  if(j > 1)
		        lines(c(rt1, rt2), c(j+.4-1, j+.95-1), lwd = mlwd, lty = 1,
                  col = "darkgrey")
			  else
		        lines(c(rt1,rt2), c(j+.4-1,j+.95-1), lwd = mlwd, lty = 1,
                  col = "darkgrey")
			}
		  }
		}
		if (plotAcrossRuns) {
 		  if (ncol(mind) <= 2) next
	      for(j in 1:(ncol(mind)-2)) {
	        ind1 <- mind[i,j]
	        ind2 <- mind[i,j+1]
	        ind3 <- mind[i,j+2]
		    if(!is.na(ind1) & is.na(ind2) & !is.na(ind3)) {
            rt1 <- rtlist[[runs[j]]][ind1]
            rt3 <- rtlist[[runs[j+2]]][ind3]
            ## cat(rt1,rtrange[1],rt3,rtrange[2],"\n")
            if( (rt1 > rtrange[1] | rt3 > rtrange[1]) &
                (rt1 < rtrange[2] | rt3 < rtrange[2]))
                lines(c(rt1,rt3), c(j+.4-1,j+.95), lwd = mlwd, lty = 1,
                      col = "green")
		    }
		  }
		  if (ncol(mind) <= 3) next
	      for(j in 1:(ncol(mind)-3)) {
	        ind1 <- mind[i,j]
	        ind2 <- mind[i,j+1]
	        ind3 <- mind[i,j+2]
	        ind4 <- mind[i,j+3]
          if(!is.na(ind1) & is.na(ind2) & is.na(ind3) & !is.na(ind4)) {
              rt1<-rtlist[[runs[j]]][ind1]
              rt4<-rtlist[[runs[j+3]]][ind4]
              ## cat(rt1,rtrange[1],rt4,rtrange[2],"\n")
              if( (rt1 > rtrange[1] | rt4 > rtrange[1]) &
                  (rt1 < rtrange[2] | rt4 < rtrange[2]))
                  lines(c(rt1,rt4), c(j+.4-1,j+1.95), lwd = mlwd, lty = 1,
                        col = "red")
		    }
		  }
	    }
	  }
  }
})





#' Add AMDIS peak detection results
#' 
#' Reads ASCII ELU-format files (output from AMDIS) and attaches them to an
#' already created \code{peaksDataset} object
#' 
#' 
#' Repeated calls to \code{parseELU} to add peak detection results to the
#' original \code{peaksDataset} object.
#' 
#' @param object a \code{peaksDataset} object.
#' @param fns character vector of same length as \code{object@rawdata} (user
#' ensures the order matches)
#' @param verbose whether to give verbose output, default \code{TRUE}
#' @param ... arguments passed on to \code{parseELU}
#' @return
#' 
#' \code{peaksDataset} object
#' @author Mark Robinson
#' @seealso \code{\link{parseELU}}, \code{\link{peaksDataset}}
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords manip
#' @examples
#' 
#' # need access to CDF (raw data) and ELU files 
#' require(gcspikelite)
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' 
#' # full paths to file names
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # create a 'peaksDataset' object and add AMDIS peaks to it
#' pd<-peaksDataset(cdfFiles[1],mz=seq(50,550),rtrange=c(7.5,8.5))
#' pd<-addAMDISPeaks(pd,eluFiles[1])
#' 
#' @export addAMDISPeaks
addAMDISPeaks<-function(object,fns=dir(,"[Eu][Ll][Uu]"),verbose=TRUE,...) {
  if(length(fns)!=length(object@rawdata))
    stop("Number of files must be the same as the number of runs (and must match).")
  mn<-1e6; mx<--1
  for(i in 1:length(object@rawrt)) {
    rng<-range(object@rawrt[[i]])
  	if(mn > rng[1]) mn<-rng[1]
	if(mx < rng[2]) mx<-rng[2]
  }
  if (verbose)
    cat("Reading retention time range:",c(mn,mx),"\n")
  for(i in 1:length(object@files)) {
    if (verbose)
	  cat("Reading",fns[i],"...")
    csin<-parseELU(fns[i],mz=object@mz,rtrange=c(mn,mx),...)
	thisind<-rep(NA,nrow(csin$tab))
	thisind.st<-rep(NA,nrow(csin$tab))
	thisind.en<-rep(NA,nrow(csin$tab))
	thisrt<-rep(NA,nrow(csin$tab))
	thisdata<-matrix(0,nrow=length(object@mz),ncol=nrow(csin$tab))
	for(j in 1:length(thisind)) {
	  thisind[j]<-which.min(abs(object@rawrt[[i]]-csin$tab[j,2]))
	  thisind.st[j]<-csin$tab[j,3]
	  thisind.en[j]<-csin$tab[j,4]
	  thisrt[j]<-object@rawrt[[i]][thisind[j]]
	  take<-which(csin$peaks[,j]>0)
	  #thisdata[take,j]<-sqrt(object@rawdata[[i]][take,thisind[j]])
	  thisdata[take,j]<-object@rawdata[[i]][take,thisind[j]]
	}
    if (verbose)
	  cat(" Done.\n")
	object@peaksind[[i]]<-thisind
	thisind.st[thisind.st==0]<-1
	object@peaksind.start[[i]]<-thisind.st
	object@peaksind.end[[i]]<-thisind.en
	object@peaksrt[[i]]<-thisrt
	object@peaksdata[[i]]<-thisdata
  }
  names(object@peaksdata)<-names(object@peaksrt)<-names(object@peaksind)<-names(object@rawdata)
  new("peaksDataset",object)
}





#' Add ChromaTOF peak detection results
#' 
#' Reads ASCII tab-delimited format files (output from ChromaTOF) and attaches
#' them to an already created \code{peaksDataset} object
#' 
#' 
#' Repeated calls to \code{parseChromaTOF} to add peak detection results to the
#' original \code{peaksDataset} object.
#' 
#' @param object a \code{peaksDataset} object.
#' @param fns character vector of same length as \code{object@rawdata} (user
#' ensures the order matches)
#' @param rtDivide number giving the amount to divide the retention times by.
#' @param verbose whether to give verbose output, default \code{TRUE}
#' @param ... arguments passed on to \code{parseChromaTOF}
#' @return
#' 
#' \code{peaksDataset} object
#' @author Mark Robinson
#' @seealso \code{\link{parseChromaTOF}}, \code{\link{peaksDataset}}
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords manip
#' @examples
#' 
#' # need access to CDF (raw data) and ChromaTOF files 
#' require(gcspikelite)
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' 
#' # full paths to file names
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' # [not run] cTofFiles<-dir(gcmsPath,"txt",full=TRUE)
#' 
#' # create a 'peaksDataset' object and add ChromaTOF peaks to it
#' pd<-peaksDataset(cdfFiles[1],mz=seq(50,550),rtrange=c(7.5,8.5))
#' # [not run] pd<-addChromTOFPeaks(pd,...)
#' 
#' @export addChromaTOFPeaks
addChromaTOFPeaks<-function(object,fns=dir(,"[Tt][Xx][Tx]"),rtDivide=60,verbose=TRUE,...) {
  if(length(fns)!=length(object@rawdata))
    stop("Number of files must be the same as the number of runs (and must match).")
  mn<-1e6; mx<--1
  for(i in 1:length(object@rawrt)) {
    rng<-range(object@rawrt[[i]])
  	if(mn > rng[1]) mn<-rng[1]
	if(mx < rng[2]) mx<-rng[2]
  }
  if (verbose)
    cat("Reading retention time range:",c(mn,mx),"\n")
  for(i in 1:length(object@files)) {
    if (verbose)
	  cat("Reading",fns[i],"...")
    csin<-parseChromaTOF(fns[i],mz=object@mz,rtrange=c(mn,mx),rtDivide=rtDivide,...)
	thisrt<-csin$tab$rt
	thisind<-rep(NA,length(thisrt))
	thisdata<-matrix(0,nrow=length(object@mz),ncol=nrow(csin$tab))
	for(j in 1:length(thisind)) {
	  thisind[j]<-which.min(abs(object@rawrt[[i]]-csin$tab$rt[j]/rtDivide))
	  thisrt[j]<-object@rawrt[[i]][thisind[j]]
	  take<-which(csin$peaks[,j]>0)
	  #thisdata[take,j]<-sqrt(object@rawdata[[i]][take,thisind[j]])
	  thisdata[take,j]<-object@rawdata[[i]][take,thisind[j]]
	}
    if (verbose)
	  cat(" Done.\n")
	object@peaksind[[i]]<-thisind
	object@peaksrt[[i]]<-thisrt
	object@peaksdata[[i]]<-thisdata
  }
  names(object@peaksdata)<-names(object@peaksrt)<-names(object@peaksind)<-names(object@rawdata)
  new("peaksDataset",object)
}



#' Plot of images of GCMS data
#' 
#' Image plots (i.e. 2D heatmaps) of raw GCMS profile data
#' 
#' For \code{peakDataset} objects, each TIC is scale to the maximum value (as
#' specified by the \code{how.near} and \code{max.near} values).  The many
#' parameters gives considerable flexibility of how the TICs can be visualized.
#' 
#' For \code{peakAlignment} objects, the similarity matrix is plotted and
#' optionally, the set of matching peaks.  \code{clusterAlignment} objects are
#' just a collection of all pairwise \code{peakAlignment} objects.
#'
#' @name plotImage
#' @aliases plotImage plotImage,peaksDataset-method
#' @param object a \code{peaksDataset} object
#' @param run index of the run to plot an image for
#' @param rtrange vector of length 2 giving start and end of the X-axis
#' (retention time)
#' @param main main title (auto-constructed if not specified)
#' @param mzrange vector of length 2 giving start and end of the Y-axis
#' (mass-to-charge ratio)
#' @param SCALE function called to scale the data (default: \code{log2})
#' @param ... further arguments passed to the \code{image} command
#' @author Mark Robinson
#' @seealso \code{\link{plot}}, \code{\link{peaksDataset}}
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords classes
#' @examples
#' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # read data
#' pd<-peaksDataset(cdfFiles[1],mz=seq(50,550),rtrange=c(7.5,8.5))
#' 
#' # image plot
#' plotImage(pd,run=1,rtrange=c(7.5,8.5),main="")
#' @importFrom gplots colorpanel
#' @importFrom graphics image
#' @export
setMethod("plotImage","peaksDataset",
          function(object, run = 1, rtrange = c(11,13), main = NULL, mzrange = c(50,200), SCALE = log2, ...){
              mz <- object@mz
  for(i in run){
    rt <- object@rawrt[[i]]
    d <- object@rawdata[[i]]	
    cols <- colorpanel(50, "green", "blue", "red")
    rind <- which(rt > rtrange[1]&rt<rtrange[2])
    mind <- which(mz>mzrange[1]&mz<mzrange[2])
    ff <- t(d[mind,rind])
    ff[ff == 0]<-min(ff[ff>0])
	if(is.null(main))
	  tit <- object@files[i]
	else
	  tit <- main
	## cat("min=",min(SCALE(ff),na.rm=TRUE)," max=",max(SCALE(ff),na.rm=T),"\n")
    image(rt[rind],mz[mind], SCALE(ff), col = cols, xlab = "retention time",
          ylab = "m/z", xlim = rtrange, main = tit, ...)
    #abline(h=c(72.5,73.5,146.5,147.5),col="white",lty=3)
    #if (!is.null(pts)) {
    #  points(pts[,1],pts[,2],pch=19)
    #  matplot(t(pts[,8:9]),t(pts[,c(2,2)]),type="l",col="black",lty=1,add=TRUE)
    #}
  }
  return(invisible(ff))
})

