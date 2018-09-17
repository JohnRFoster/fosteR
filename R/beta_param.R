##---------------------------------------------------------------------------------##
##              Beta Distribution Parameter estimation: Method of Moments          ##
##---------------------------------------------------------------------------------##

#' Beta Distribution
#'
#' This function takes a vector with data that is [0, 1]
#' Uses mean and variance of input vector to estimate
#' the shape parameters (alpha & beta) for the beta distribution
#' @param vec vector with values between 0 and 1
#' @export
#' @examples beta_param(vec)

beta_param <- function(vec){
  mu <- mean(vec)
  var <- var(vec)
  mom <- ((mu*(1-mu))/var) - 1
  a <- mu*mom
  b <- (1-mu)*mom
  param <- data.frame(alpha = a, beta = b)
  return(param)
}
