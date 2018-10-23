##---------------------------------------------##
##             Logit Transformation            ##
##---------------------------------------------##

#' Logit Transformation
#'
#' This function takes a number and transforms it to
#' the logit scale
#'
#' @param p Number that you wish to transform
#' @export
#' @examples logit(p)


logit <- function(p){
  param <- 1/(1 + exp(-p))
  return(param)
}
