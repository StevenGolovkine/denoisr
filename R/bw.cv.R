######################################################################################
#                    Functions for bandwiths cross-validation                        #
######################################################################################
library(np)
library(tidyverse)

#' Perform bandwith selection by cross-validation for the kernel regression estimator.
#' 
#' @param data List of curves to estimate by kernel regression
#' @return A vector containing the mimimum, the median, the mean and the maximum 
#' bandwith selected by cross-validation for every curve.
bw.cv <- function(data){
  # Get the number of curves
  N <- length(data)
  
  # Get a sample of the curves of size (at least) log(N)
  sample_data <- sample(data, ceiling(log(N)))
  
  # Perform a kernel regression estimator for each of the sample data
  bw_list <- sample_data %>% 
    map_dbl(~ sqrt(5) * np::npregbw(.x$x ~ .x$t, 
                         regtype = 'll', # Local Linear Regression
                         bwmethod = 'cv.ls', # Least Square Cross Validation
                         ckertype = 'epanechnikov')$bw # Kernel
    )
  
  bw <- c(min(bw_list), median(bw_list), mean(bw_list), max(bw_list))
  return(bw)
}

