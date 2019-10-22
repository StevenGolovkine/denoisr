###############################################################################
#        Functions for L0 parameter estimation using regularity               #
###############################################################################
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
estimate.L0 <- function(data, sigma = NULL, H0 = 0, 
                        k0 = 2, t0 = 0, density = NULL){
  
  S_N <- data
  
  # Estimate mu
  mu_hat <- S_N %>% map_int(~ length(.x$t)) %>% mean()
  
  # Estimate L0
  theta <- function(v, k, idx) (v[idx + 2*k - 1] - v[idx + k])**2
  eta <- function(v, k, idx, H) (v[idx + 2*k - 1] - v[idx + k])**(2*H)
  
  nume <- 0; deno <- 1
  if (is.null(density)) {# Case where the density is not known
    if (is.null(sigma)) {# Subcase where sigma is not known
      idxs <- S_N %>%
        map_dbl(~ min(order(abs(.x$t - t0))[seq_len(4*k0 - 2)]))
      a <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = 2*k0 - 1, idx = .y)) %>% 
        mean()
      b <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>% 
        mean()
      c <- S_N %>% 
        map2_dbl(idxs, ~ eta(.x$t, k = 2*k0 - 1, idx = .y, H = H0)) %>% 
        mean()
      d <- S_N %>% 
        map2_dbl(idxs, ~ eta(.x$t, k = k0, idx = .y, H = H0)) %>% 
        mean()
      if ((a - b > 0) & (c - d > 0)) {
        nume <- a - b
        deno <- c - d
      }
    } else {# Subcase where sigma is known
      idxs <- S_N %>% 
        map_dbl(~ min(order(abs(.x$t - t0))[seq_len(2*k0)]))
      a <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>% 
        mean()
      b <- S_N %>% 
        map2_dbl(idxs, ~ eta(.x$t, k = k0, idx = .y, H = H0)) %>% 
        mean()
      if ((a - 2*sigma**2 > 0) & b > 0) {
        nume <- a - 2*sigma**2
        deno <- b
      }
    }
  } else {# Case where the density is known (only the uniform case)
    if (is.null(sigma)) {# Subcase where sigma is not known
      idxs <- S_N %>%
        map_dbl(~ min(order(abs(.x$t - t0))[seq_len(4*k0 - 2)]))
      a <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = 2*k0 - 1, idx = .y)) %>% 
        mean()
      b <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>% 
        mean()
      if (a - b > 0) nume <- a - b
      deno <- (2 ** (2 * H0) - 1) * ((k0 - 1) / (mu_hat + 1))**(2 * H0) 
    } else {# Subcase where sigma is known
      idxs <- S_N %>% 
        map_dbl(~ min(order(abs(.x$t - t0))[seq_len(2*k0)]))
      a <- S_N %>% 
        map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>% 
        mean()
      if (a - 2*sigma**2 > 0) nume <- a - 2*sigma**2
      deno <- ((k0 - 1) / (mu_hat + 1))**(2 * H0)
    }
  }
  
  L0_hat <- (nume / deno)**0.5 
  
  return(L0_hat)
}
