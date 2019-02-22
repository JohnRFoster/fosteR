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
#' @param site Site. One of "Green Control","Henry Control","Tea Control"
#' @param data Data used in model run.
#' @param params matrix of parameters from mcmc: as.matrix(out$params)
#' @param ic matrix of states estiamted from mcmc: as.matrix(out$predict)
#' @param Nmc Number of random draws from mcmc output
#' @param draw Vector of randomly sampled row numbers from mcmc output
#' @export
#' @examples predict_state_one("parameter", "Green Control", as.matrix(out$params),as.matrix(out$predict),500,2,sample.int(nrow(as.matrix(out$params)), 500, replace = TRUE))


predict_state_one <- function(type, site, params, ic, data, Nmc, draw){

  # call data from model run for storage dimensions and indexing
  if(site == "Green Control"){
    N_est <- data$N_est[1]
    N_days <- data$N_days[1]
    dt.index <- data$dt.index[1,]
    df <- data$df[1,]
  } else if(site == "Henry Control"){
    N_est <- data$N_est[2]
    N_days <- data$N_days[2]
    dt.index <- data$dt.index[2,]
    df <- data$df[2,]
  } else {
    N_est <- data$N_est[3]
    N_days <- data$N_days[3]
    dt.index <- data$dt.index[3,]
    df <- data$df[3,]
  }

  # storage
  pred <- array(dim = c(3,N_est,Nmc))
  A <- array(0, dim=c(3,3,N_days))

  # mean parameters
  param.mean <- apply(params, 2, mean)

  # select appropriate initial conditions
  if("ic" %in% type){
    IC <- ic[draw,]
  } else {
    IC.mean <- apply(ic, 2, mean)
    IC.names <- names(IC.mean)
    IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
    colnames(IC) <- IC.names
  }

  # select appropriate parameters
  if("parameter" %in% type){
    phi.l.mu <- params[draw,"phi.l.mu"]
    phi.n.mu <- params[draw,"phi.n.mu"]
    phi.a.mu <- params[draw,"phi.a.mu"]
    grow.ln.mu <- params[draw,"grow.ln.mu"]
    grow.na.mu <- params[draw,"grow.na.mu"]
    repro.mu <- params[draw,"repro.mu"]
  } else {
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
  }

  # select appropraite covariance matrix
  if("process" %in% type){
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

    # convert from precision to standard dev
    for(i in 1:Nmc){
      SIGMA[,,i] <- solve(SIGMA[,,i])
    }
  } else {
    SIGMA <- array(0, dim = c(3,3,Nmc))
  }

  # run mcmc sampling
  for(m in 1:Nmc){

    # draw transition matrix
    A[1,1,] <- logit(phi.l.mu[m])
    A[1,3,] <- exp(repro.mu[m])
    A[2,1,] <- logit(grow.ln.mu[m])
    A[2,2,] <- logit(phi.n.mu[m])
    A[3,2,] <- logit(grow.na.mu[m])
    A[3,3,] <- logit(phi.a.mu[m])

    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){                 # loop over the number of sampling days - 1
      if(t == 1){
        TRANS <- A[,,1] %*% A[,,2]
        for(d in 3:df[1]){
          TRANS <- TRANS %*% A[,,d]
        }
      } else {
        for(d in 1:df[t]){
          if(d == 1){
            TRANS <- A[,,dt.index[t-1]] %*% A[,,dt.index[t-1]+1]
          } else {
            TRANS <- TRANS %*% A[,,dt.index[t-1]+d]
          }
        }
      }

      ## initial condition
      l <- IC[m, paste("x[1,",t,"]",sep="")]
      n <- IC[m, paste("x[2,",t,"]",sep="")]
      a <- IC[m, paste("x[3,",t,"]",sep="")]
      obs <- as.matrix(c(l,n,a),3,1)
      Ex <- TRANS %*% obs
      est.mvnorm <- rmvnorm(1,Ex,SIGMA[,,m])
      pred[1,t,m] <- max(est.mvnorm[1], 0)
      pred[2,t,m] <- max(est.mvnorm[2], 0)
      pred[3,t,m] <- max(est.mvnorm[3], 0)
    }
  }
  return(pred)
} # close function
