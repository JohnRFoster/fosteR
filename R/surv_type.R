##---------------------------------------------------------------------------------##
##                            Tick Survival Data intake                            ##
##---------------------------------------------------------------------------------##

#' Tick Survival Data Intake function
#'
#' This funciton reads the raw csv file containing tick survival data
#' from the Cary institute survival experiments (soil cores) and transforms
#' the data to be fed to JAGS
#'
#' @param file Path to CSV
#' @param subset one of "Fed Larvae", "Flat Larvae", "Fed Nymph", "Flat Nymph", "Overwintering Nymph"
#' @param season default = NULL, must set to "winter or "summer for larvae
#' @export
#' @examples surv.type(file, "Fed Nymph")


surv_type <- function(file, subset, season = NULL){

  rawdata <- read.csv(file)
  dat <- rawdata
  dat$Bin <- 1:nrow(dat) # create column that is just row numbers

  dat$Deployment_Name <- as.character(dat$Deployment_Name)
  dat$Site <- as.numeric(dat$Site)  # convert site to numeric IDs

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

  if(dat$Deployment_Name == "Flat Larvae" ||
     dat$Deployment_Name == "Flat Nymph" ||
     dat$Deployment_Name == "Overwintering Nymph"){

    # JAGS data for flat ticks
    data = list(y = cbind(dat$N_Deployed, dat$N_Recovered),
                N_days = N_Days,
                site.index = dat$Site,
                N = nrow(dat),
                sites = unique(dat$Site),
                N_site = length(unique(dat$Site)))
  } else {

    # JAGS data for fed ticks (includes number successfully molted)
    data = list(y = cbind(dat$N_Deployed, dat$N_Survived),
                N_days = N_Days,
                site.index = dat$Site,
                N = nrow(dat),
                sites = unique(dat$Site),
                N_site = length(unique(dat$Site)))
  }
  return(data)
}




