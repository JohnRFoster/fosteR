##-------------------------------------------------------------------##
##                 Split Tick States from JAGS output                ##
##-------------------------------------------------------------------##

#' Split state estimates from JAGS Daily model into respective sites
#'
#' This function takes the state matrix from the HB tick model and splits the sites
#' @param state Matric of state estimates containing multiple sites
#' @param n.sites The number of sites to split
#' @examples split_tick_states(state, 3)
#' @export


split_tick_states <- function(state,n.sites){
  IC.list <- list()
  seq <- 1:n.sites
  for(i in 1:n.sites){
    IC <- state[,grep(paste(i, "]", sep = ""), colnames(state))]
    IC.list[[i]] <- IC
  }
  return(IC.list)
}
