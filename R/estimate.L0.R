######################################################################################
#               Functions for L0 parameter estimation using regularity               #
######################################################################################
library(tidyverse)

#' Perform the estimation of L0. 
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param H0 An estimation of H0.
#' @param k The sampling time to consider (>= 2).
#' @param density Density of the sampling points (currently, only consider uniform 
#'        sampling points).
#' @param sigma True value of sigma
#'         If null, change estimate
#' @return An estimation of L0.
estimate.L0 <- function(data, H0 = 0, k = 2, density = NULL, sigma=NULL){
  
  S_N <- data
  
  # Estimate mu
  mu_hat <- S_N %>% map_int(~ length(.x$t)) %>% mean()
  
  # Estimate L0
  substract_double_element <- function(v, i){
    return((v[4 * i - 3] - v[2 * i - 1])**2 - (v[2 * i - 1] - v[i])**2)
  }
  substract_double_time <- function(v, i, H0){
    return(abs(v[4 * i - 3] - v[2 * i - 1])**(2 * H0) - abs(v[2 * i - 1] - v[i])**(2 * H0))
  }
  substract_first_element <- function(v, i) (v[2 * i - 1] - v[i])**2
  substract_first_time <- function(v, i, H0) abs(v[2 * i - 1] - v[i])**(2 * H0)
  substract_sigma <- function(nb, sigma=0) nb - 2*sigma**2
  
  if(is.null(density)){
    if(is.null(sigma)){
      nume <- S_N %>%
        map_dbl(~ substract_double_element(.x$x, k)) %>%
        mean()
      deno <- S_N %>%
        map_dbl(~ substract_double_time(.x$t, k, H0)) %>%
        mean()
    } else{
      nume <- S_N %>%
        map_dbl(~ substract_first_element(.x$x, k)) %>%
        mean() %>%
        substract_sigma(sigma)
      deno <- S_N %>%
        map_dbl(~ substract_first_time(.x$t, k, H0)) %>%
        mean()
    }
  } else{
    if(is.null(sigma)){
      nume <- S_N %>%
        map_dbl(~ substract_double_element(.x$x, k)) %>%
        mean()
      deno <- (2 ** (2 * H0) - 1) * ((k - 1) / (mu_hat + 1))**(2 * H0) 
    } else{
      nume <- S_N %>%
        map_dbl(~ substract_first_element(.x$x, k)) %>%
        mean() %>%
        substract_sigma(sigma)
      deno <- (2 ** (2 * H0) - 1) * ((k - 1) / (mu_hat + 1))**(2 * H0)
    }
  }
  
  L0_hat <- (nume / deno)**0.5 
  
  return(L0_hat)
}
