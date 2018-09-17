##---------------------------------------------------------------------------------##
##                           ANNUAL Cary Met Data intake                           ##
##---------------------------------------------------------------------------------##

#' Annual mean monthly RH, Temp, and Precip
#'
#' Function to create Annual Met from Cary weather stations
#'
#' @param start start year
#' @param end end year
#' @param path file path to csv if not in working directory
#' @param center Should that data be centered? default is TRUE
#' @export
#' @examples met_caryAnnual(1995, 2005)

met_caryAnnual <- function(start, end, path = "", center = TRUE){

  file <- "Met_Cary.csv"
  met <- read.csv(paste(path, file, sep = ""))
  met <- met[, c("DATE","MAX_TEMP","MAX_RH","TOT_PREC")]           # select daily max temp, max rh, precip
  met.y <- as.data.frame(strsplit(as.character(met$DATE),"-"))     # split date by "-"
  met.y <- t(met.y[1,])                                            # transpose date and select year
  met.y <- apply(met.y, 2, as.integer)                             # make integer
  met <- cbind(met.y, met[,-1])                                    # cbind year to met data
  colnames(met) <- c("year","MAX_TEMP","MAX_RH","TOT_PREC")
  year <- seq(start, end)                                                # create vector of years
  annual.met <- data.frame()
  for(i in 1:length(year)){                                        # loop over years and calculate mean
    xx <- met[met$year==year[i],]
    for(c in 1:ncol(met)){
      annual.met[i, c] <- mean(xx[,c], na.rm = TRUE)
    }
  }
  colnames(annual.met) <- colnames(met)
  if(center == FALSE){                                             # if center = FALSE, return data
    return(annual.met)
  }
  if(center == TRUE){
    annual.met <- apply(annual.met[,-1], 2, scale, scale = FALSE)    # center data around mean
    annual.met <- cbind(year, annual.met)                            # cbind year
    return(annual.met)
  }
}

