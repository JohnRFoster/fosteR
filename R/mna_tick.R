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

mna_tick_current <- function(sites = c("Green Control","Henry Control","Tea Control")){
  minpositive <- function(x) min(x[x > 0])
  mna <- list()
  for(i in 1:length(sites)){
    tick <- tick_cary(sites[i], "individual")
    ch <- suppressWarnings(ch_cary(sites[i]))
    ks <- known_states(ch)
    tick.dates <- as.Date(colnames(tick))
    tick.seq <- seq.Date(tick.dates[1],tick.dates[length(tick.dates)],1)
    mouse.dates <- as.Date(colnames(ks))
    mna.obs <- apply(ks, 2, sum)
    mna.x <- vector()
    for(j in 1:length(tick.seq)){
      close.date <- which(tick.seq[j]-mouse.dates == minpositive(tick.seq[j]-mouse.dates))
      mna.x[j] <- mna.obs[close.date]
    }
    mna[[i]] <- mna.x
  }
  out <- matrix(NA, 3, length(mna[[which.max(lengths(mna))]]))
  for(i in 1:length(sites)){
    for(j in 1:length(mna.test[[i]]))
      out[i,j] <- mna.test[[i]][j]
  }
  return(out)
}
