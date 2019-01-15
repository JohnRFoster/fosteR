##---------------------------------------------------------------------------------##
##                          One step ahead prediction - null                       ##
##---------------------------------------------------------------------------------##

#' Predict next time step from current - null model (demographic params + random effects)
#'
#' This funciton takes output from the null mcmc tick model and makes predictions one
#' step at a time. It also partitions uncertainty into initial condition ("ic"), parameter
#' ("parameter"), process ("process"), and random effects ("random effect")
#'
#' @param type Type of uncertainty to partition, one of "deterministic", "ic", "parameter", "process", or "random effects"
#' @param params matrix of parameters from mcmc: as.matrix(out$params)
#' @param ic matrix of states estiamted from mcmc: as.matrix(out$predict)
#' @param Nmc Number of random draws from mcmc output
#' @param s Site. Numeric 1-3. Maps to "Green Control","Henry Control","Tea Control"
#' @param draw Vector of randomly sampled row numbers from mcmc output
#' @param sites Vector of sites used in model run, default are control sites: c("Green Control","Henry Control","Tea Control")
#' @export
#' @examples predict_state_one("parameter",as.matrix(out$params),as.matrix(out$predict),500,2,sample.int(nrow(as.matrix(out$params)), 500, replace = TRUE))


predict_state_one <- function(type, params, ic, Nmc, s, draw,
                          sites = c("Green Control","Henry Control","Tea Control")){

  # call data from model run for storage dimensions and indexing
  data <- cary_ticks_JAGS(sites)
  N_est <- data$N_est
  df <- data$df

  # storage
  pred <- array(dim = c(3,N_est[s],Nmc))
  A <- matrix(0,3,3)

  # mean parameters
  param.mean <- apply(params, 2, mean)

  # select appropriate initial conditions
  if(type == "deterministic"){
    IC <- t(as.matrix(apply(ic, 2, mean)))
  } else {
    IC <- ic[draw,]
  }

  # select appropriate parameters
  if(type == "deterministic" | type == "ic"){
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
  } else {
    phi.l.mu <- params[draw,"phi.l.mu"]
    phi.n.mu <- params[draw,"phi.n.mu"]
    phi.a.mu <- params[draw,"phi.a.mu"]
    grow.ln.mu <- params[draw,"grow.ln.mu"]
    grow.na.mu <- params[draw,"grow.na.mu"]
    repro.mu <- params[draw,"repro.mu"]
  }

  # select appropriate random effect error
  if(type == "random effect"){
    tau.11 <- 1/sqrt(params[draw,"tau.11"])
    tau.22 <- 1/sqrt(params[draw,"tau.22"])
    tau.33 <- 1/sqrt(params[draw,"tau.33"])
    tau.13 <- 1/sqrt(params[draw,"tau.13"])
    tau.21 <- 1/sqrt(params[draw,"tau.21"])
    tau.32 <- 1/sqrt(params[draw,"tau.32"])
    alpha.11 <- rnorm(Nmc, 0, tau.11)
    alpha.22 <- rnorm(Nmc, 0, tau.22)
    alpha.33 <- rnorm(Nmc, 0, tau.33)
    alpha.13 <- rnorm(Nmc, 0, tau.13)
    alpha.21 <- rnorm(Nmc, 0, tau.21)
    alpha.32 <- rnorm(Nmc, 0, tau.32)
  } else if(type == "process" | type == "parameter") {
    alpha.11 <- params[draw,paste("alpha.11[",s,"]",sep = "")]
    alpha.22 <- params[draw,paste("alpha.22[",s,"]",sep = "")]
    alpha.33 <- params[draw,paste("alpha.33[",s,"]",sep = "")]
    alpha.13 <- params[draw,paste("alpha.13[",s,"]",sep = "")]
    alpha.21 <- params[draw,paste("alpha.21[",s,"]",sep = "")]
    alpha.32 <- params[draw,paste("alpha.32[",s,"]",sep = "")]
  } else {
    alpha.11 <- rep(param.mean[paste("alpha.11[",s,"]",sep = "")], Nmc)
    alpha.22 <- rep(param.mean[paste("alpha.22[",s,"]",sep = "")], Nmc)
    alpha.33 <- rep(param.mean[paste("alpha.33[",s,"]",sep = "")], Nmc)
    alpha.13 <- rep(param.mean[paste("alpha.13[",s,"]",sep = "")], Nmc)
    alpha.21 <- rep(param.mean[paste("alpha.21[",s,"]",sep = "")], Nmc)
    alpha.32 <- rep(param.mean[paste("alpha.32[",s,"]",sep = "")], Nmc)
  }

  # select appropraite covariance matrix
  if(type == "deterministic" | type == "ic" | type == "parameter"){
    SIGMA <- array(0, dim = c(3,3,Nmc))
  } else {
    SIGMA <- array(NA, dim = c(3,3,Nmc))
    SIGMA[1,1,] <- params[draw,"SIGMA[1,1]"]
    SIGMA[1,2,] <- params[draw,"SIGMA[1,2]"]
    SIGMA[1,3,] <- params[draw,"SIGMA[1,3]"]
    SIGMA[2,1,] <- params[draw,"SIGMA[2,1]"]
    SIGMA[2,2,] <- params[draw,"SIGMA[2,2]"]
    SIGMA[2,3,] <- params[draw,"SIGMA[2,3]"]
    SIGMA[3,1,] <- params[draw,"SIGMA[3,1]"]
    SIGMA[3,2,] <- params[draw,"SIGMA[3,2]"]
    SIGMA[3,3,] <- params[draw,"SIGMA[3,3]"]
  }

  # run mcmc sampling
  for(m in 1:Nmc){

    # draw transition matrix
    A[1,1] <- logit(phi.l.mu[m] + alpha.11[m])
    A[1,3] <- exp(repro.mu[m] + alpha.13[m])
    A[2,1] <- logit(grow.ln.mu[m] + alpha.21[m])
    A[2,2] <- logit(phi.n.mu[m] + alpha.22[m])
    A[3,2] <- logit(grow.na.mu[m] + alpha.32[m])
    A[3,3] <- logit(phi.a.mu[m] + alpha.33[m])

    ## aggrigate transition matricies
    for(t in 1:(N_est[s]-1)){
      TRANS <- A %*% A
      for(d in 2:df[s,t]){
        TRANS <- TRANS %*% A
      }

      ## initial condition
      l <- IC[m, paste("x[1,",t,",",s,"]",sep="")]
      n <- IC[m, paste("x[2,",t,",",s,"]",sep="")]
      a <- IC[m, paste("x[3,",t,",",s,"]",sep="")]
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      pred[1:3,t,m] <- rmvnorm(1,Ex,solve(SIGMA[,,m]))
    }
  }
  return(pred)
} # close function
