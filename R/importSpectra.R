##' Read the mass spectra from an external msp file
##'
##' Read the mass spectra from an external file in msp format. The format is
##' used in NIST search library database.
##' @title importSpec
##' @param file a .msp file from NIST search library database
##' @return list conaining the mass spctra
##' @author riccardo.romoli@unifi.it
##' @export
importSpec <- function(file){
    ## read msp lib
    lib <- scan(file, what = "", sep = "\n", quiet = TRUE)
    ## separate each mass spec
    starts <- which(regexpr("[Nn][Aa][Mm][Ee]:", lib) == 1)
    ends <- c(starts[-1] - 1, length(lib))
    ## loop to extract the mass spec into a list
    list.spec <- lapply(1:length(starts), function(z){
        ## meta data
#        browser()
        comp <- lib[starts[z]:ends[z]]
        numPeaks.idx <- which(regexpr("[Nn][Uu][Mm] [Pp][Ee][Aa][Kk][Ss]:", comp) == 1)
        metaData <- comp[1:numPeaks.idx - 1]
        md <- strsplit(metaData, split = ": ")
        md1 <- sapply(md, "[[", 1)
        md2 <- sapply(md, "[", 2)
        metaData.list <- setNames(as.list(md2), md1)
        ## mass spec
        nlines <- length(comp)
        npeaks <- as.numeric(strsplit(comp[numPeaks.idx], ":")[[1]][2])
        peaks.idx <- (numPeaks.idx + 1):nlines
        pks <- gsub("^ +", "", unlist(strsplit(comp[peaks.idx], ";")))
        ## il separatore potrebbe essere anche (), va trovata una soluzione
        pks.2 <- gsub("^ +", "", unlist(strsplit(comp[peaks.idx], "[()]")))
gsub("\"", "", pks.2)

        pks <- pks[pks != ""]
        if (length(pks) != npeaks)
        {
            stop("Not the right number of peaks in compound", metaData.list$NAME)
        }
        ## error due to the presence of tab \t instead of a blank space
        if(length(grep(pattern = "\t", pks)) > 0)
        {
            pklst <- strsplit(pks, "\t")
            pklst <- lapply(pklst, function(x) x[x != ""])
        }
        else if(length(grep(pattern = "\\s", pks)) > 0)
        {
            pklst <- strsplit(pks, " ")
            pklst <- lapply(pklst, function(x) x[x != ""])
        }
        else
        {
            cat("Some formatting errors in the msp library \n")
        }

        mz <- as.numeric(sapply(pklst, "[[", 1))
        mz <- round(mz)
        int <- as.numeric(sapply(pklst, "[[", 2)) # error
        ##
        finaltab <- matrix(c(mz, int), ncol = 2)

        if (any(table(mz) > 1))
        {
            warning("Duplicate mass in compound ", metaData.list$NAME,
                    " (CAS ", metaData.list$NAME, ")... summing up intensities")
            finaltab <- aggregate(finaltab[,2],
                                  by = list(finaltab[,1]),
                                  FUN = sum)
        }
        colnames(finaltab) <- c("mz", "intensity")
        c(metaData.list, list(spec = finaltab))
    }
    )
    return(list.spec)
}


##' Calculate the distance between a reference mass spectrum
##'
##' Calculate the distance between a reference mass spectrum and one from the
##' sample
##' @title matchSpec
##' @param spec1 reference mass spectrum
##' @param outList the return of \code{\link{gatherInfo}}
##' @param whichSpec the entry number of outList
##' @return the distance between the reference mass spectrum and the others
##' @author Riccardo Romoli
##' @export
matchSpec <- function(spec1, outList, whichSpec){
  ## if(whichSpec == 143) browser()
  ## first get the average spec from the gatherList
  averageInt <-
    apply(outList[[whichSpec]]$data, MARGIN = 1, FUN = mean, na.rm = TRUE)
  ## outList[[whichSpec]]$data[,1:49] # si ma perch49?!?
  ## normalize the intensity
  normInt <- sapply(averageInt, function(x){
    i <- which.max(averageInt)
    (100*x)/averageInt[i]
  }
  )
  ## combine mz and normalized intensity
  pspec <- matrix(c(outList[[whichSpec]]$mz, normInt), ncol = 2)
  colnames(pspec) <- c("mz", "intensity")
  libSpec <- spec1$spec
  ## merge the spectra to the same mz range and remove the NA
  mztomerge <- data.frame(mz = min(c(min(libSpec[,"mz"], na.rm = TRUE),
                                     min(pspec[,"mz"], na.rm = TRUE))) :
                            max(c(max(libSpec[,"mz"], na.rm = TRUE),
                                  max(pspec[,"mz"], na.rm = TRUE)))
                          )
  pspec.merged <- merge(pspec, mztomerge, by = "mz", all = TRUE)
  libSpec.merged <- merge(libSpec, mztomerge, by = "mz", all = TRUE)
  pspec.merged[is.na(pspec.merged)] <- 0
  libSpec.merged[is.na(libSpec.merged)] <- 0
  ## calculate the distance among the spectra
  distance <- normDotProduct(as.matrix(pspec.merged[,"intensity"]),
                             as.matrix(libSpec.merged[,"intensity"]))
  return(distance)
}


##' The function calculate the distance between each mas spec in the msp file
##' and the aligned mass spec from each sampe
##'
##' Return the distance matrix
##' @title distToLib
##' @param mspLib a .msp file from NIST
##' @param outList an object from gatherInfo()
##' @return the distance matrix between the mass spec and the aligned spec
##' @author Riccardo Romoli
##' @export
distToLib <- function(mspLib, outList ){
  mspDist <- lapply(1:length(mspLib),
                    function(x)
                    {
                      ## trace
                      cat(x, "\n")
                      sapply(1:length(outList),
                             function(y)
                             {
                               ## cat("y is = ", y, "\n") # for debug only
                               matchSpec(mspLib[[x]], outList,
                                         whichSpec = y)
                             }
                             )
                    })
  mtx.dist <- do.call(rbind, mspDist)
  rownames(mtx.dist) <- sapply(mspLib, "[[", 1)
  colnames(mtx.dist) <- sapply(1:length(outList),
                               function(x)
                               {
                                 paste0("outListFold", x)
                               }
                               )
  return(mtx.dist)
}


##' The head-to-tail-plot for the mass spectra
##'
##' Head-to-tail-plot to visually compare the mass spectra
##' @title Head to tail plot
##' @param specFromLib the mass spectra obtained from the .msp file
##' @param specFromList the mass spectra obtained from \code{\link{gatherInfo}}
##' @return the plot
##' @author Riccardo Romoli
##' @export
headToTailPlot <- function(specFromLib, specFromList){
    libSpec <- specFromLib$spec
    pspec <- specFromList
    ## get average and normalized intensity of the mass spec
    averageInt <- apply(pspec$data, MARGIN = 1, FUN = mean, na.rm = TRUE)
    ## normalize the intensity
    normInt <- sapply(averageInt, function(x){
        i <- which.max(averageInt)
        (1000*x)/averageInt[i]
    }
    )
    pspec.av <- data.frame(mz = pspec$mz, intensity = normInt)
    ## merge the spectra to the same mz range and remove the NA
    mztomerge <- data.frame(mz = min(c(min(libSpec[,"mz"], na.rm = TRUE),
                                       min(pspec.av[,"mz"], na.rm = TRUE))) :
                              max(c(max(libSpec[,"mz"], na.rm = TRUE),
                                    max(pspec.av[,"mz"], na.rm = TRUE)))
                            )
    pspec.merged <- merge(pspec.av, mztomerge, by = "mz", all = TRUE)
    libSpec.merged <- merge(libSpec, mztomerge, by = "mz", all = TRUE)
    pspec.merged[is.na(pspec.merged)] <- 0
    libSpec.merged[is.na(libSpec.merged)] <- 0
    ## now the plot
    plot(pspec.merged, type = "h", ylim = c(-1000, 1000), main = "Head to Tail Plot")
    points(libSpec.merged[,"mz"], y = -libSpec.merged[,"intensity"]*10, col = 2,
           type = "h")
}
