##---------------------------------------------------------------------------------##
##                            Tick Survival Data intake                            ##
##---------------------------------------------------------------------------------##

#' Tick Survival Data Intake function
#'
#' This funciton reads the raw csv file containing tick survival data
#' from the Cary institute survival experiments (soil cores)
#' @param file Path to CSV
#' @param season "winter" or "summer"
#' @param winter.start date when winter season starts "yyyy-mm-dd"
#' @param winter.end date when winter season ends "yyyy-mm-dd"
#' @param life.stage "Larvae" or "Nymph"
#' @param flat.fed "Flat" or "Fed"
#' @export
#' @examples surv.type(file, "summer", "2017-10-01", "2018-05-31", "Nymph", "Flat")


surv_type <- function(file, season, winter.start, winter.end, life.stage, flat.fed){

  rawdata <- read.csv(file)
  dat <- rawdata

  dat$Site <- as.numeric(dat$Site)  # convert site to numeric IDs

  # convert from factor to date
  dat$Date_Deployed <- as.Date(as.character(dat$Date_Deployed), "%m/%d/%Y")
  dat$Date_Retrieved <- as.Date(as.character(dat$Date_Retrieved), "%m/%d/%Y")

  start.date <- dat$Date_Deployed   # extract start date
  end.date <- dat$Date_Retrieved    # extract end date

  # rows the represent over-winter experiments
  over.winter <- which(end.date <= as.Date(winter.end) &
                         start.date >= as.Date(winter.start))

  # rows the represent over-summer experiments
  over.summer <- which(end.date < as.Date(winter.start) |
                         start.date > as.Date(winter.end))

  ## over winter or within-season experiments?
  if(season == "winter"){
    dat <- dat[over.winter,]   # all over winter data
  }
  if(season == "summer"){
    dat <- dat[over.summer,]   # all summer data
  }

  # dat <- subset(dat, N_Recovered > 0)      # subset only cores where ticks recovered
  dat <- subset(dat, Stage == life.stage)       # subset life stage
  dat <- subset(dat, Flat_Fed == flat.fed)      # subset flat or fed ticks
  dat <- subset(dat, Animal_Disturbance == "N") # subset no animal disutbance
  dat <- subset(dat, Tick_Source == "NY")       # subset only ticks from NY

  dat <- dat[,c("Site",
                "Date_Deployed",
                "Date_Retrieved",
                "N_Deployed",
                "N_Recovered",
                "N_Recovered_Engorged",
                "N_Recovered_Molted",
                "N_Molted_After",
                "N_Survived")]

  dat[is.na(dat)] <- 0  # convert remaining NAs to 0

  N_Molted_Total <- vector()
  for(j in 1:nrow(dat)){
    N_Molted_Total[j] <- dat$N_Recovered_Molted[j] + dat$N_Molted_After[j]
  }

  dat <- cbind(dat, N_Molted_Total)

  # Number of days in soil
  N_Days <- difftime(dat$Date_Retrieved, dat$Date_Deployed)

  dat <- cbind(dat, N_Days)

  return(dat)
}
