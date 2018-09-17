##---------------------------------------------------------------------------------##
##                    Tick life-stage CIs for an individual site                   ##
##---------------------------------------------------------------------------------##


#' CI's for tick state estimates by life stage from an HB model (includes all sites)
#' from Cary Institute
#'
#' @param state matrix of iterations from bayes sampling
#' @param site What site to create CI's? One of "Green", "Henry", or "Tea"
#' @param nsamp How many random samples to draw? Default = 5000
#' @export
#' @examples ci_TickStatesHB(jags.out, "Tea")

ci_TickStatesHB <- function(state, site, nsamp = 5000){

  ## this function extracts state iterations from a
  ## specified site, takes nsamp random samples
  ## and creates CI for the three life stages,
  ## then returns one object (life stage CI for
  ## an individual site across time)

  if(site == "Green"){
    ss <- state[,grep("1]", colnames(state))]
  }
  if(site == "Henry"){
    ss <- state[,grep("2]", colnames(state))]
  }
  if(site == "Tea"){
    ss <- state[,grep("3]", colnames(state))]
  }

  end <- ncol(ss)

  l.seq <- seq(1,end,by = 3) # sequence for larvae columns
  n.seq <- seq(2,end,by = 3) # sequence for nymph columns
  a.seq <- seq(3,end,by = 3) # sequence for adult columns

  larv <- ss[,l.seq]         # matrix with larve iterations
  nymph <- ss[,n.seq]        # matrix with nymph iterations
  adult <- ss[,a.seq]        # matrix with adult iterations

  samp <- sample.int(nrow(larv),nsamp) # random samples

  # ci for the three life stages
  ci.larvae <- apply(larv[samp,],2,quantile,c(0.025,0.5,0.975))
  ci.nymph <- apply(nymph[samp,],2,quantile,c(0.025,0.5,0.975))
  ci.adult <- apply(adult[samp,],2,quantile,c(0.025,0.5,0.975))

  ci <- rbind(ci.larvae, ci.nymph, ci.adult)     # combine to one object

  rownames(ci) <- c("larvae 2.5%", "larvae 50%", "larvae 97.5%",
                    "nymph 2.5%", "nymph 50%", "nymph 97.5%",
                    "adult 2.5%", "adult 50%", "adult 97.5%")
  return(ci)
}
