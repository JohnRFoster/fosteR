##---------------------------------------------------------------------------------##
##                          Monthly Cary Met Data intake                           ##
##---------------------------------------------------------------------------------##

#' Monthly mean monthly RH, Temp, and Precip
#'
#' Function to create Monthly Met from Cary weather stations
#'
#' @param start start month (yyyymm)
#' @param end end month (yyyymm)
#' @param path file path to csv if not in working directory
#' @param center Should that data be centered? default is TRUE
#' @export
#' @examples met_caryAnnual(199505, 200511)

met_caryMonthly <- function(start, end, path = "", center = TRUE){
  file <- "Met_Cary.csv"
  met <- read.csv(paste(path, file, sep = ""))
  met <- met[, c("DATE","MAX_TEMP","MAX_RH","TOT_PREC")]           # select daily max temp, max rh, precip
  met.y <- as.data.frame(strsplit(as.character(met$DATE),"-"))
  met.y <- t(met.y[-3,])
  met.y <- paste(met.y[,1], met.y[, 2], sep = "")
  met <- cbind(met.y, met[,-1])
  colnames(met) <- c("yr_mon","MAX_TEMP","MAX_RH","TOT_PREC")
  yr_mon <- unique(met.y)
  yr_mon.index <- yr_mon[which(yr_mon == start):which(yr_mon == end)]
  monthly.met <- data.frame()
  for(i in 1:length(yr_mon.index)){
    xx <- met[met$yr_mon == yr_mon.index[i], -1]
    xx <- apply(xx, 2, mean, na.rm = TRUE)
    monthly.met <- rbind(monthly.met, xx)
  }
  monthly.met <- cbind(yr_mon.index, monthly.met)
  colnames(monthly.met) <- colnames(met)
  if(center == FALSE){                                             # if center = FALSE, return data
    return(monthly.met)
  }
  if(center == TRUE){
    monthly.met <- apply(monthly.met[,-1], 2, scale, scale = FALSE)    # center data around mean
    monthly.met <- cbind(yr_mon.index, monthly.met)                            # cbind year
    monthly.met <- apply(monthly.met, 2, as.numeric)
    return(monthly.met)
  }
}
##----------------------##
##     End Function     ##
##----------------------##
