################################################################################
#       Functions for sigma parameter estimation using regularity              #
################################################################################


#' Perform the estimation of sigma 
#' 
#' @param data List of curves to estimate by kernel regression.
#' @return An estimation of sigma.
estimate.sigma <- function(data){
  
  S_N <- data
  
  # Estimate sigma
  sigma_hat <- S_N %>% estimateSigma()
  
  return(sigma_hat)
}
