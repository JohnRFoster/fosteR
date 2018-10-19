##---------------------------------------------------------------------------------##
##              Tick Survival Data intake - AA ticks from Brian Allan              ##
##---------------------------------------------------------------------------------##


#' Tick Survival Data Intake function for Brian Allan's AA ticks
#'
#' This funciton reads the raw csv file containing tick survival data
#' from the BA's survival experiments (soil cores) and transforms
#' the data to be fed to JAGS
#'
#' @param life.stage "Nymph" or "Adult"
#' @param path File path to csv if not in working directory
#' @export
#' @examples Allan_surv("Nymph")


Allan_surv <- function(path = "", life.stage){

  if(life.stage == "Nymph"){

  nymph <- read.csv(paste(path,"Tick_Mortality_Nymph.csv",
                          sep = ""), header = FALSE)
  nymph <- nymph[,-c(1:4)]  # extract just tick numbers
  date <- apply(nymph[1,], 2, as.character) # convert first row to character
  date <- as.Date(date, "%m/%d/%Y") # date of census
  nymph <- apply(nymph[-1,], 2, as.numeric) # data
  N_days <- length(seq.Date(date[1], date[length(date)], 1)) # total number of days
  days.diff <- as.numeric(diff(date)) # days between census's

  data <- list(y = nymph,
               N = nrow(nymph),
               census = ncol(nymph),
               N_days = N_days,
               delta.days = days.diff)

  return(data)


  } else if(life.stage == "Adult"){

  adult <- read.csv(paste(path,"Tick_Mortality_Adult.csv",
                          sep = ""),header = FALSE)
  adult <- adult[,-c(1:4)]  # extract just tick numbers
  date <- apply(adult[1,], 2, as.character) # convert first row to character
  date <- as.Date(date, "%m/%d/%Y") # date of census
  adult <- apply(adult[-1,], 2, as.numeric) # data
  adult <- adult[-6,] # remove site 6; data inconsistent
  N_days <- length(seq.Date(date[1], date[length(date)], 1)) # total number of days
  days.diff <- as.numeric(diff(date)) # days between census's

  data <- list(y = adult,
               N = nrow(adult),
               census = ncol(adult),
               N_days = N_days,
               delta.days = days.diff)

  return(data)
  }
}



