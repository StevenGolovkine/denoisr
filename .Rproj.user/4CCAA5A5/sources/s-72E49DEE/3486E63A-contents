######################################################################################
#               Functions for H0 parameter estimation using regularity               #
######################################################################################
library(tidyverse)


#' Perform the estimation of H0. 
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param sigma True value of sigma
#'         If null, change estimate
#' @param k0 For the computation of the gap between the different observations. Should
#'         be set as k0 = M / max(8, log(M)).
#' @param t0 The starting time for the estimation of H0. We consider the 8k - 7 nearest
#'         points of t0 for the estimation of H0 when sigma is unknown.
#' 
#' @return An estimation of H0.
estimate.H0 <- function(data, sigma = NULL, k0 = 2, t0 = 0){
  
  S_N <- data
  
  # Estimate H_0
  theta <- function(v, t, k, idx) (v[idx + 2*k - 1] - v[idx + k])**2
  substract_sigma <- function(nb, sigma=0) nb - 2*sigma**2
  
  two_log_two <- 2*log(2)
  if(is.null(sigma)){ # Case where sigma is unknown
    idxs <- S_N %>%
      map_dbl(~ min(order(abs(.x$t - t0))[seq_len(8*k0 - 6)]))
    a <- S_N %>% 
      map2_dbl(idxs, ~ theta(.x$x, .x$t, k = 4*k0 - 3, idx = .y)) %>% 
      mean()
    b <- S_N %>% 
      map2_dbl(idxs, ~ theta(.x$x, .x$t, k = 2*k0 - 1, idx = .y)) %>% 
      mean()
    c <- S_N %>% 
      map2_dbl(idxs, ~ theta(.x$x, .x$t, k = k0, idx = .y)) %>%
      mean()
    first_part <- log(a - b)
    second_part <- log(b - c)
    
  } else{ # Case where sigma is known
    idxs <- S_N %>%
      map_dbl(~ min(order(abs(.x$t - t0))[seq_len(2*k0 + 1)]))
    first_part <- S_N %>% 
      map2_dbl(idxs, ~ theta(.x$x, .x$t, k = k0 + 1, idx = .y)) %>% 
      mean() %>% 
      substract_sigma(sigma) %>%
      max(c(., 1e-5)) %>%
      log()
    second_part <- S_N %>% 
      map2_dbl(idxs, ~ theta(.x$x, .x$t, k = k0, idx = .y)) %>% 
      mean() %>% 
      substract_sigma(sigma) %>% 
      max(c(., 1e-5)) %>%
      log()
  }

  H0_hat <- (first_part - second_part) / two_log_two

  return(H0_hat)
}
