##---------------------------------------------------------------------------------##
##                        Mouse MNA closest to tick sampling                       ##
##---------------------------------------------------------------------------------##

#' Mouse MNA closest to tick sampling
#'
#' Function that finds the minimum number alive of mice nearest in time to
#' tick sampling date, returns a vector of MNA
#' @param tick.date vector of dates from tick sampling
#' @param mouse.date vecor of dates from mouse captures
#' @param known.states known states matrix
#' @export
#' @example mna_tick()

mna_tick <- function(tick.date, mouse.date, known.states){
  # tick.date = vector of tick sampling dates
  # mouse.date = vector of mice sampling dates
  # known.states = known states matrix

  # this function matches the mna for a given site
  # to the closest tick sampling date in the past

  mna <- apply(known.states, 2, sum)
  tick.mna <- vector()
  for(t in 1:length(dates.tick)){
    xx <- which(abs(dates.tick[t]-dates.mouse) == min(abs(dates.tick[t]-dates.mouse)))
    if(dates.mouse[xx] > dates.tick[t]){
      yy <- dates.mouse[1:(xx-1)]
      zz <- which(abs(dates.tick[t]-yy) == min(abs(dates.tick[t]-yy)))
      aa <- as.character(dates.mouse[zz])
    } else {
      aa <- as.character(dates.mouse[xx])
    }
    tick.mna[t] <- mna[aa][[1]]
  }
  return(tick.mna)
}
