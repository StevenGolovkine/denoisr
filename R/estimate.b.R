######################################################################################
#           Functions for bandwith parameter estimation using regularity             #
######################################################################################
library(tidyverse)


#' Perform the estimation of the bandwith 
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param sigma An estimation of sigma.
#' @param H0 An estimation of H0.
#' @param L0 An estimation of L0.
#' @param kernel Which kernel to use?
#' @return An estimation of H0.
estimate.b <- function(data, sigma=0, H0=0.5, L0=1, K='epanechnikov'){
  
  S_N <- data
  
  # Set kernel constants
  if (K == 'epanechnikov') {
    K_norm2 <- 0.6
    K_phi <- 3 / ((H0 + 1) * (H0 + 3))
  } else if (K == 'beta') {
    K_norm2_f <- function(x, alpha = 1, beta = 1){x**(2*(alpha - 1)) * (1 - x)**(2*(beta - 1)) / beta(alpha, beta)**2}
    phi <- function(x, H0, alpha = 1, beta = 1){abs(x**(alpha - 1) * (1 - x)**(beta - 1)) * x**H0 / abs(beta(alpha, beta))}
    K_norm2 <- integrate(K_norm2_f, lower = lower, upper = upper)$value
    K_phi <- integrate(phi, lower = 0, upper = 1, H0 = H0)$value
  } else{
    K_norm2 <- 1
    K_phi <- 1 / (H0 + 1)
  }
  
  # Estimate mu
  mu_hat <- S_N %>% map_int(~ length(.x$t)) %>% mean()
  
  # Estimate b
  #H0_tilde <- H0 - (1 / log(mu_hat))
  nume <- (sigma**2 * K_norm2 * factorial(floor(H0)))
  deno <- (H0 * L0 * K_phi)
  frac <- nume / deno
  b_hat <- (frac / mu_hat)**(1 / (2*H0 + 1))
  
  return(b_hat)
}

#' Perform the estimation of the bandwith using CV
#' 
#' @param data list of curves to estimate by kernel regression
estimate.b.cv <- function(data){
  require(np); require(doParallel); require(dplyr)
  
  cl <- parallel::detectCores() %>%  - 1 %>% parallel::makeCluster()
  doParallel::registerDoParallel(cl)
  
  bw_list <- foreach(j = 1:length(data)) %dopar% {
    sqrt(5) * np::npregbw(x ~ t, data = data[[j]], 
                          bwmethod = 'cv.ls', # Least Square Cross Validation
                          ckertype = 'epanechnikov', # Kernel used
                          regtype = 'lc')$bw # Local Constant Regression
  }
  
  parallel::stopCluster(cl)
  
  return(mean(unlist(bw_list)))
}
