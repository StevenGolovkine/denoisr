######################################################################################
#                Functions for bandwiths estimation using regularity                 #
######################################################################################
library(tidyverse)

Rcpp::sourceCpp('./src/estimate_sigma.cpp')

#' Perform bandwith selection by estimation of regularity of the curves.
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param pct Percentage of data to consider for the estimation.
#' @param sigma The value of the standard deviation of the noise.
#'               If sigma=NULL, it will be estimated.
#' @param l The sampling time to consider (> 1).
#' @param K Kernel used (default: 'epanechnkov')
#' 
#' @return An estimation of the bandwith
estimate.bw <- function(data, pct=0.1, sigma=NULL, s_0=NULL, L_0=NULL, l=2, K='epanechnikov', 
                        est_list=c('mu' = TRUE, 'sigma' = TRUE, 's_0' = TRUE, 'L_0' = TRUE, 'b' = TRUE, 'risk' = TRUE)){
  
  # Step 0 - Set kernel constante
  if(K == 'epanechnikov'){
    K_norm2 <- 0.6
    phi <- function(x, s_0) {0.75 * (1 - x**2) * abs(x)**s_0}
    upper <- 1
    lower <- -1
  } else if(K == 'beta'){
    K_norm2_f <- function(x, alpha = 1, beta = 1){x**(2*(alpha - 1)) * (1 - x)**(2*(beta - 1)) / beta(alpha, beta)**2}
    phi <- function(x, s_0, alpha = 1, beta = 1){abs(x**(alpha - 1) * (1 - x)**(beta - 1)) * x**s_0 / abs(beta(alpha, beta))}
    upper <- 1
    lower <- 0
  } else{
    K_norm2 <- 1
    phi <- function(x, s_0) {1}
    upper <- 1
    lower <- 0
  }
  
  # Step 1 - Split the dataset into 2 parts.
  N_0 <- round(length(data)*pct, 0)
  idx <- sample(1:length(data), N_0, replace = TRUE)
  S_N_0 <- data[idx]
  S_N <- data[-idx]
  
  # Step 2 - Estimate mu
  mu_hat <- 0
  if(est_list['mu']){
    mu_hat <- S_N_0 %>% map_int(~ length(.x$t)) %>% mean()
  }

  # Step 2 (bis) - Estimate sigma
  sigma_hat <- sigma
  if(est_list['sigma']){
    if(is.null(sigma)) sigma_hat <- S_N_0 %>% estimateSigma()
  }

  # Step 3 - Estimate s_0
  substract_first_element <- function(v, i) (v[i] - v[1])**2
  substract_sigma <- function(nb, sigma=0) nb - 2*sigma**2
  
  # We assume that the first element of each curve is at time 0.
  s_0_hat <- s_0
  if(est_list['s_0']){
    if(is.null(s_0)){
      two_log_two <- 2*log(2)
      first_part <- S_N_0 %>% 
        map_dbl(~ substract_first_element(.x$x, 3)) %>% 
        mean() %>% 
        substract_sigma(sigma_hat) %>% 
        #abs() %>% # Take the absolute value to get rid of NaN value.
        log()
      second_part <- S_N_0 %>% 
        map_dbl(~ substract_first_element(.x$x, 2)) %>% 
        mean() %>% 
        substract_sigma(sigma_hat) %>% 
        #abs() %>% # Take the absolute value to get rid of NaN value.
        log()
      s_0_hat <- (first_part - second_part) / two_log_two
    }
  }

  # Step 4 - Estimate L_0
  L_0_hat <- L_0
  if(est_list['L_0']){
    if(is.null(L_0)){
      nume <- S_N_0 %>%
        map_dbl(~ substract_first_element(.x$x, l)) %>%
        mean() %>%
        substract_sigma(sigma_hat) %>%
        #abs() # Take the absolute value to get rid of NaN value.
      deno <- ((l - 1) / (mu_hat + 1))**(2 * s_0_hat)
      L_0_hat <- (nume / deno)**0.5 
    }
  }
  
  # Step 5 - Estimate b
  b_hat <- 0
  if(est_list['b']){
    if(K == 'beta'){
      K_norm2 <- integrate(K_norm2_f, lower = lower, upper = upper)$value
    }
    # Step 5(i) - Computation of the constante
    s_0_tilde <- s_0_hat #- (1 / log(mu_hat))
    nume <- (sigma_hat**2 * K_norm2 * factorial(floor(s_0_tilde)))
    deno <- (s_0_tilde * L_0_hat * integrate(phi, lower = lower, upper = upper, s_0 = s_0_tilde)$value)
    c <- nume / deno 
    
    # Step 5(ii) - Computation of the bandwith
    b_hat <- (c / mu_hat)**(1 / (2*s_0_hat + 1))
  }
  
  # Step 6 - Computation of the risk on S_n
  risk <- list()
  if(est_list['risk']){
    risk <- S_N %>% estimateRisk(b = b_hat)
  }
  
  return(list(mu_hat = mu_hat,
              sigma_hat = sigma_hat,
              s_0_hat = s_0_hat,
              L_0_hat = L_0_hat,
              b_hat = b_hat,
              risk = risk))
}
