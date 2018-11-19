##---------------------------------------------------------------------------------##
##                    Tick Data intake; 3 x n sampling matrix                      ##
##---------------------------------------------------------------------------------##

#' Tick Data Intake
#'
#' Raw tick data was not great, so the following function just takes the cleaned data
#' and subsets for a specified grid and data type
#' @param path File path to csv, default = ""
#' @param grid What grid do you want data? One of Green Control, Henry Control, Tea Control, Green Experimental, Henry Experimental, Tea Experimental
#' @param type How should the population counted? One of density, individual, or delta
#' @export
#' @examples tick_cary(grid = "Green Control", type = "individual")

tick_cary <- function(grid, type, path = ""){

  file <- "tick_cleaned"
  t <- read.csv(paste(path, file, sep = ""))                   # read file
  t <- subset(t, Grid == grid)                                # subset grid
  d <- t$DATE
  if(type == "individual"){                                     # subset the number of ticks caught
    t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
  }
  if(type == "density"){
    t <- t[,c("Larvae.m2", "Nymphs.m2", "Adults.m2")]
  }
  if(type == "delta"){
    t <- t[,c( "delta.larvae","delta.nymph", "delta.adult")]
  }
  if(type != "individual" && type != "density" && type != "delta") {
    print("ERROR: Did you forget to set type? Must be one of density, individual, or delta")
  }
  dat.t <- t(t)                                               # transpose
  dat.t <- apply(dat.t, 2, as.integer)
  colnames(dat.t) <- d
  return(dat.t)
}

