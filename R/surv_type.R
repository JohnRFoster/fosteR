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
  dat$Bin <- 1:nrow(dat) # create column that is just row numbers

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

  #dat <- subset(dat, N_Recovered > 0)      # subset only cores where ticks recovered
  dat <- subset(dat, Stage == life.stage)       # subset life stage
  dat <- subset(dat, Flat_Fed == flat.fed)      # subset flat or fed ticks
  dat <- subset(dat, Animal_Disturbance == "N") # subset no animal disutbance
  dat <- subset(dat, Tick_Source == "NY")       # subset only ticks from NY

  # sequence of all days within range of subsetted data
  days.seq <- seq.Date(min(dat$Date_Deployed),
                       max(dat$Date_Retrieved), 1)

  # Number of days in soil
  N_Days <- as.numeric(difftime(dat$Date_Retrieved, dat$Date_Deployed))

  # create matrix
  core.mat <- matrix(0,nrow(dat),length(days.seq))
  colnames(core.mat) <- as.character(days.seq)

  # convert date to character
  deploy.date <- as.character(dat$Date_Deployed)
  retrieve.date <- as.character(dat$Date_Retrieved)

  # fill in core matrix and create indexing variable for start date
  start.index <- vector(length = length(deploy.date))
  for(t in 1:length(deploy.date)){
    n.deploy <- dat[t,"N_Deployed"]                  # number deployed
    deploy.col <- which(deploy.date[t]==days.seq)    # date in sequence it was deployed
    core.mat[t,deploy.col] <- n.deploy               # assign to corresponding column

    n.recovered <- dat[t,"N_Recovered"]                # number recovered
    recovered.col <- which(retrieve.date[t]==days.seq) # date in seq it was recovered
    core.mat[t,recovered.col] <- n.recovered           # assing to corresponding column

    start.index[t] <- deploy.col   # indexing variable for deploy date
  }

  # JAGS data for flat ticks
  if(flat.fed == "Flat"){
    data = list(y = core.mat,
                start.index = start.index,
                N_days = N_Days,
                site.id = dat$Site,
                N = nrow(dat))
  }
  # JAGS data for fed ticks (includes number successfully molted)
  if(flat.fed == "Fed"){
    data = list(y = core.mat,
                start.index = start.index,
                N_days = N_Days,
                site.id = dat$Site,
                N = nrow(dat),
                N_Molted = dat$N_Survived)
  }

  return(data)
}



