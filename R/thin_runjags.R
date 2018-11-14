##---------------------------------------------------------------------------------##
##                                  Thin run.jags                                  ##
##---------------------------------------------------------------------------------##

#' Thin run.jags output to specified interval
#'
#' This function reads in the output (.RData) file containing a run.jags object
#' and thins the mcmc chains to the user specified interval. Currently supports
#' 3 chains only!
#'
#' @param RData File path to .RData object containg chains
#' @param thin Total number of samples to keep
#' @export
#' @examples thin_runjags(RData, thin)

thin_runjags <- function(RData, thin){
  load(RData)
  start <- as.integer(rownames(run.jags[["mcmc"]][[1]])[1])
  end <- dim(run.jags[["mcmc"]][[1]])[1]
  seq <- seq(start, end, length = thin)
  mcmc <- list()
  mcmc[[1]] <- run.jags[["mcmc"]][[1]][seq,]
  mcmc[[2]] <- run.jags[["mcmc"]][[2]][seq,]
  mcmc[[3]] <- run.jags[["mcmc"]][[3]][seq,]
  return(mcmc)
}
