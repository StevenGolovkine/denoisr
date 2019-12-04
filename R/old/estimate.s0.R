#' ######################################################################################
#' #               Functions for H0 parameter estimation using regularity               #
#' ######################################################################################
#' library(tidyverse)
#' 
#' 
#' #' Perform the estimation of H0. 
#' #' 
#' #' @param data List of curves to estimate by kernel regression.
#' #' @param sigma True value of sigma
#' #'         If null, change estimate
#' #' @param t0 The starting time for the estimation of H0.
#' #' 
#' #' @return An estimation of H0.
#' estimate.s0 <- function(data, sigma=NULL, t0=0, method = 'forward'){
#'   
#'   S_N <- data
#'   
#'   # Estimate s_0
#'   theta <- function(v, t, k, t0, method =  'forward'){
#'     idx <- sum(t < t0)
#'     if(method == 'backward') idx <- idx - (2 * k + 1)
#'     return((v[2*k - 1] - v[k])**2)
#'   } 
#'   substract_element <- function(v, t, i, t0, method){
#'     idx <- sum(t <= t0)
#'     if(method == 'backward') return((v[idx - 1] - v[idx -1 - 2*i])**2 - (v[idx - 1] - v[idx - 1 - i])**2)
#'     else return((v[2*i + 1 + idx] - v[1 + idx])**2 - (v[i + 1 + idx] - v[1 + idx])**2)
#'   } 
#'   substract_first_element <- function(v, t, i, t0){
#'     idx <- sum(t <= t0)
#'     return((v[i + idx] - v[1 + idx])**2)
#'   }
#'   substract_sigma <- function(nb, sigma=0) nb - 2*sigma**2
#'   
#'   two_log_two <- 2*log(2)
#'   if(is.null(sigma)){
#'     first_part <- S_N %>% 
#'       map_dbl(~ theta(.x$x, .x$t, k = 5, t0 = t0, method = method) - theta(.x$x, .x$t, k = 3, t0 = t0, method = method)) %>% 
#'       #map_dbl(~ substract_element(.x$x, .x$t, 2, t0 = t0, method = method)) %>% 
#'       mean() %>% 
#'       max(c(., 1e-5)) %>%
#'       log()
#'     second_part <- S_N %>% 
#'       map_dbl(~ theta(.x$x, .x$t, k = 3, t0 = t0, method = method) - theta(.x$x, .x$t, k = 2, t0 = t0, method = method)) %>%
#'       #map_dbl(~ substract_element(.x$x, .x$t, 1, t0 = t0, method = method)) %>% 
#'       mean() %>% 
#'       max(c(., 1e-5)) %>%
#'       log()
#'   } else{
#'     first_part <- S_N %>% 
#'       map_dbl(~ theta(.x$x, .x$t, k = 3, t0 = t0, method = method)) %>% 
#'       mean() %>% 
#'       substract_sigma(sigma) %>%
#'       max(c(., 1e-5)) %>%
#'       log()
#'     second_part <- S_N %>% 
#'       map_dbl(~ theta(.x$x, .x$t, k = 2, t0 = t0, method = method)) %>% 
#'       mean() %>% 
#'       substract_sigma(sigma) %>% 
#'       max(c(., 1e-5)) %>%
#'       log()
#'   }
#' 
#'   H0_hat <- (first_part - second_part) / two_log_two
#' 
#'   return(H0_hat)
#' }
