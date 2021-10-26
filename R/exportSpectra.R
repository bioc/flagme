#' exportSpectra
#'
#' Write the mass spectum into a .msp file to be used in NIST search.
#'
#' Write the mass spectum into a .msp file to be used in NIST search.
#'
#' @param object an object of class "peaksDataset"
#' @param outList an object created using the gatherInfo() function
#' @param spectra numeric. The number of the mass spectra to be printed. It
#' correspond to the number of the peak in the plot() and the number of the
#' peak in the gatherInfo() list.
#' @param normalize logical. If the mass spectra has to be normalized to 100
#' @return a .msp file
#' @author riccardo.romoli@@unifi.com
#' @export exportSpectra
exportSpectra <- function (object, outList, spectra, normalize = TRUE)
{
    spectra <- as.numeric(spectra)
    mz <- outList[[spectra]]$mz
    abu <- outList[[spectra]]$data
    if(normalize)
    {
        for(i in 1:ncol(abu)){
            if(is.na(sum(abu[, i])))
                next
            abu[, i] <- 100 * abu[, i]/abu[which.max(abu[, i]), i]
        }
    }
    ## specnum <- table(apply(abu, 2, sum))
    ## now I have a number of spectrum that is equal to the number of time that
    ## the mass spectum was match accross the samples. As a first approach
    ## I decide to print only the first mass spec.
    idx <- 1
    spec <- paste(mz, abu[,idx])
    msp <- rbind(paste("NAME: Variable", spectra),
                 paste("COMMENT:", round(object@peaksrt[[idx]][spectra], digits = 2), "min"),
                 paste("FORMULA:"), paste("MW:"), paste("CAS:"), paste("SYNONYM:"),
                 paste("Num Peaks:", length(spec)), matrix(unlist(spec), nrow = length(unlist(spec)), ncol = 1))

    write(msp, file = paste0(spectra,".msp"), sep = "\n")

}
