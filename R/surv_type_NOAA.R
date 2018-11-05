##---------------------------------------------------------------------------------##
##                            Tick Survival Data intake                            ##
##---------------------------------------------------------------------------------##

#' Tick Survival Data Intake function
#'
#' This funciton is a wrapper around the surv_type function that also reads in
#' NOAA data and reads the raw csv file containing tick survival data
#' from the Cary institute survival experiments (soil cores) and transforms
#' the data to be fed to JAGS
#'
#' @param file Path to CSV
#' @param subset one of "Fed Larvae", "Flat Larvae", "Fed Nymph", "Flat Nymph", "Overwintering Nymph"
#' @param season default = NULL, must set to "winter or "summer for larvae
#' @export
#' @examples surv.type(file, "Fed Nymph")


surv_type_NOAA <- function(file, subset, season = NULL){

  rawdata <- read.csv(file)
  dat <- rawdata
  dat$Bin <- 1:nrow(dat) # create column that is just row numbers

  dat$Deployment_Name <- as.character(dat$Deployment_Name)

  # convert from factor to date
  dat$Date_Deployed <- as.Date(as.character(dat$Date_Deployed), "%m/%d/%Y")
  dat$Date_Retrieved <- as.Date(as.character(dat$Date_Retrieved), "%m/%d/%Y")

  start.date <- dat$Date_Deployed   # extract start date
  end.date <- dat$Date_Retrieved    # extract end date

  dat <- subset(dat, Animal_Disturbance == "N") # subset no animal disutbance
  dat <- subset(dat, Tick_Source == "NY")       # subset only ticks from NY
  dat <- subset(dat, !is.na(dat$N_Recovered))   # subset ticks with conflicting info in N_Recovered
  dat <- subset(dat, Deployment_Name == subset) # subset data to use

  # parse larve into winter and summer subset
  if(dat$Deployment_Name == "Flat Larvae" || dat$Deployment_Name == "Fed Larvae"){
    if(season == "winter"){
      dat <- subset(dat, Date_Retrieved >= "2018-01-01")
    } else {
      dat <- subset(dat, Date_Retrieved < "2018-01-01")
    }
  }

  # sequence of all days within range of subsetted data
  days.seq <- seq.Date(min(dat$Date_Deployed),
                       max(dat$Date_Retrieved), 1)

  # Number of days in soil
  N_Days <- as.numeric(difftime(dat$Date_Retrieved, dat$Date_Deployed))

  # read in NOAA data
  met.cl <- read.csv("CampLejune_NOAA.csv")
  met.fd <- read.csv("FortDrum_NOAA.csv")
  met.wp <- read.csv("WestPoint_NOAA.csv")

  # set rownames to date
  row.names(met.cl) <- met.cl[,2]
  row.names(met.fd) <- met.fd[,2]
  row.names(met.wp) <- met.wp[,2]

  # only keep met variables
  met.cl <- met.cl[,3:5]
  met.fd <- met.fd[,3:5]
  met.wp <- met.wp[,3:5]

  met.cl$tmax <- scale(met.cl$tmax, scale = FALSE)
  met.cl$tmin <- scale(met.cl$tmin, scale = FALSE)
  met.fd$tmax <- scale(met.fd$tmax, scale = FALSE)
  met.fd$tmin <- scale(met.fd$tmin, scale = FALSE)
  met.wp$tmax <- scale(met.wp$tmax, scale = FALSE)
  met.wp$tmin <- scale(met.wp$tmin, scale = FALSE)

  # all met array
    # dim 1 = date
    # dim 2 = met variable
    # dim 3 = site; "CL","FD","WP"
  met.all <- array(data = c(unlist(met.cl),unlist(met.fd),unlist(met.wp)),
                   dim = c(nrow(met.cl),ncol(met.cl),3),
                   dimnames = list(rownames(met.cl),
                                   colnames(met.cl),
                                   c("CL", "FD", "WP")))

  met.x <- array(data = NA,dim = c(3,nrow(dat),length(days.seq)),
                 dimnames = list(c("Tmax","Tmin","precip"),
                                 as.character(dat$TC_ID),
                                 as.character(days.seq)))

  ## met.x dimensions:
    # dim 1 = met variable; index = m
    # dim 2 = soil cores; index = c
    # dim 3 = date; index = t
  for(c in 1:nrow(dat)){
    seq <- seq.Date(dat$Date_Deployed[c],dat$Date_Retrieved[c], 1)
    start <- which(seq[1] == days.seq) # start index for column
    end <- which(seq[length(seq)] == days.seq) # end index for column
    for(m in 1:ncol(met.cl)){
      for(t in start:end){
        site.dim <- which(dat$Site[c] == dimnames(met.all)[[3]])
        row.dim <- which(days.seq[t] == dimnames(met.all)[[1]])
        met.x[m,c,t] <- met.all[row.dim,m,site.dim]
      }
    }
  }

  # indexing variable for when in sequence core was deployed
  start.index <- vector()
  for(i in 1:nrow(dat)){
    start.index[i] <- which(dat$Date_Deployed[i] == days.seq)
  }

  dat$Site <- as.numeric(dat$Site)  # convert site to numeric IDs

  if(dat$Deployment_Name == "Flat Larvae" ||
     dat$Deployment_Name == "Flat Nymph" ||
     dat$Deployment_Name == "Overwintering Nymph"){

    # JAGS data for flat ticks
    data = list(y = cbind(dat$N_Deployed, dat$N_Recovered),
                N_days = N_Days,
                site.index = dat$Site,
                N = nrow(dat),
                start.index = start.index,
                sites = unique(dat$Site),
                N_site = length(unique(dat$Site)),
                met = met.x)
  } else {

    # JAGS data for fed ticks (includes number successfully molted)
    data = list(y = cbind(dat$N_Deployed, dat$N_Survived),
                N_days = N_Days,
                site.index = dat$Site,
                N = nrow(dat),
                start.index = start.index,
                sites = unique(dat$Site),
                N_site = length(unique(dat$Site)),
                met = met.x)
  }
  return(data)
}


