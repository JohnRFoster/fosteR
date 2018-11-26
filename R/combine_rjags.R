##---------------------------------------------------------------------------------##
##                                  Combine rjags                                  ##
##---------------------------------------------------------------------------------##

#' Combine separate and segemented rjags chains
#'
#' This function reads segemented rjags chains that are saved in sequential .rds files,
#' combines them to form a single chain, and then combines multiple chains to give an
#' mcmc object for model diagnostics.
#'
#' @param sites character vector of sites to estimate
#' @export
#' @examples


combine_rjags <- function(path, chain.num, start = 1, end, thin = NULL){
  n.chain <- length(chain.num)
  chains <- list()
  for(c in 1:n.chain){
    for(i in start:end){
      iter <- i-start+1
      file <- paste(path, chain.num[c], sep = "")
      mcmc[iter] <- read.RDS(file = paste(file, i, "rds", sep = "."))
    }
    chains[[c]] <- combine.mcmc(mcmc)
  }
  if(!is.null(thin)){
    seq <- round(seq(1, by = dim(chains[[1]][1]/thin,length = thin)))
    for(i in 1:n.chain){
      chains[[c]] <- as.mcmc(chains[[c]][seq,])
    }
  }
  jags.out <- as.mcmc(chains)

}


