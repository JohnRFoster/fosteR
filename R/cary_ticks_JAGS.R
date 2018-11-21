##---------------------------------------------------------------------------------##
##                           Daily Cary tick data for JAGS                         ##
##---------------------------------------------------------------------------------##

#' Daily Cary tick data for JAGS
#'
#' This function reads the raw tick drag data from the Cary institute and formats it to be
#' ready for JAGS. This data is for running models to estimate the state (abundance of
#' ticks) for each life stage every seven days and for each sampling occasion for the
#' entire time series
#'
#' @param sites character vector of sites to estimate
#' @param state.interval Default = NULL, the interval to estimate the state (7 = weekly)
#' @export
#' @examples cary_ticks_JAGS(sites = c("Green Control","Henry Control","Tea Control"))


cary_ticks_JAGS <- function(sites, state.interval = NULL){
  N_site <- length(sites)               # number of sites
  raw.dat <- read.csv("tick_cleaned")   # read in data
  raw.dat$DATE <- as.Date(raw.dat$DATE) # convert to date

  # storage
  date.seq <- list(length(sites)) # all days within sampling window
  N_days <- vector()              # total number of days within sampling window
  df <- list()                    # number of days between sampling events
  dt.index <- list()              # number of days elaplsed from the first sampling occasion
  tick <- list(length(sites))     # data
  samp.dates <- list()            # sampling days
  interval.seq <- list()          # days to estimate latent state
  N_est <- vector()               # total number of days to estimate the latent state

  # separate by site
  for(j in 1:length(sites)){
    t <- subset(raw.dat, Grid == sites[j])
    t <- t[,c("DATE", "n_larvae", "n_nymphs", "n_adults")]
    tick[[j]] <- t
    samp.dates[[j]] <- t$DATE
    date.seq[[j]] <- data.frame(seq.Date(t$DATE[1], t$DATE[nrow(t)], 1))
    N_days[j] <- nrow(date.seq[[j]])

    # if we want to estimate the state between sampling occasions
    if(!is.null(state.interval)){
      interval <- seq(1,by=state.interval,length=floor(nrow(date.seq[[j]]) / state.interval))
      interval.seq[[j]] <- sort(as.Date(c(samp.dates[[j]],date.seq[[j]][interval,]),format = "%Y-%m-%d"))
      interval.seq[[j]] <- unique(interval.seq[[j]])
      N_est[j] <- length(interval.seq[[j]])
      df[[j]] <- as.numeric(diff(interval.seq[[j]]))
      index <- c(df[[j]],0)
      dt.index[[j]] <- cumsum(index)

    # if we want to estimate the state on just the sampling days
    } else {
      N_est[j] <- nrow(t)
      interval.seq[[j]] <- samp.dates[[j]]
      df[[j]] <- as.numeric(diff(t$DATE))
      index <- c(df[[j]],0)
      dt.index[[j]] <- cumsum(index)
    }
  }

  df.mat <- matrix(NA,length(sites),max(sapply(df,length)))
  for(j in 1:length(sites)){
    for(t in 1:length(df[[j]]))
    df.mat[j,t] <- df[[j]][[t]]
  }
  dt.index.mat <- matrix(NA,length(sites),max(sapply(dt.index,length)))
  for(j in 1:length(sites)){
    for(t in 1:length(dt.index[[j]]))
      dt.index.mat[j,t] <- dt.index[[j]][[t]]
  }

  # all tick array
    # dim 1 = life stage (larvae, nymph, adult)
    # dim 2 = date
    # dim 3 = site
  all.tick <- array(NA,dim = c(3,max(N_est),N_site))

  # match tick date to interval.seq, place in array
  for(s in 1:N_site){
    for(i in 1:nrow(tick[[s]])){
      larv <- tick[[s]][i,"n_larvae"]
      nymph <- tick[[s]][i,"n_nymphs"]
      adult <- tick[[s]][i,"n_adults"]
      col.index <- which(tick[[s]]$DATE[i] == interval.seq[[s]])
      all.tick[1,col.index,s] <- larv
      all.tick[2,col.index,s] <- nymph
      all.tick[3,col.index,s] <- adult
    }
  }

  data <- list(y = all.tick,
               dt.index = dt.index.mat,
               df = df.mat,
               N_est = N_est,
               N_days = N_days)

 return(data)

}
