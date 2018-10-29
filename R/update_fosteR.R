##---------------------------------------------------------------------------------##
##                                Automatic update                                 ##
##---------------------------------------------------------------------------------##

#' Automatic update function
#'
#' @export
#' @examples update_fosteR()



update_fosteR <- function(){
  library(devtools)
  install_github("JohnRFoster/fosteR")
}
