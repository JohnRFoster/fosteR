##---------------------------------------------------------------------------------##
##                                Automatic update                                 ##
##---------------------------------------------------------------------------------##

#' Automatic update function
#'
#' @export
#' @examples update_fosteR()



update_fosteR <- function(force = FALSE){
  library(devtools)
  install_github("JohnRFoster/fosteR",force = force)
}
