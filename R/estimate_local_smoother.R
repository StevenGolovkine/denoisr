################################################################################
#         Functions for covariance estimation with local linear smoother       #
################################################################################

#' Estimation of the mean with local smoother
#' 
#' This function performs the estimation of the mean function based on the paper
#' by Yao, Müller and Wang.
#' 
#' @param data A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param time A vector containing the times where the mean is computed.
#' @param regtype A character string specifying which type of kernel regression
#'  estimator to use. lc specifies a local-constant estimator (Nadaraya-Watson) 
#'  and ll specifies a local-linear estimator. Defaults to lc.
#' 
#' @return A vector representing the mean curve.
#'  
#' @references Yao, Müller and Wang - Functional Data Analysis for Sparse
#'  Longitudinal Data (2005) - Journal of the American Statistical Association
#' @export
local_mean <- function(data, time, regtype = 'lc'){
  data_ <- list2cai(data)
  
  band <- np::npregbw(x ~ time, data = data_,
                  ckertype = 'epanechnikov', regtype = regtype)
  mod <- np::npreg(band, newdata = data.frame(time = time))
  mod$mean
}

#' Estimation of the covariance with local smoother
#' 
#' This function performs the estimation of the covariance function based on the
#' paper by Cai and Yuan.
#' 
#' @param data A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param t Vector, contains the times where the covariance is computed.
#' @param s Vector, contains the times where the covariance is computed on the
#'  other axis. IF NULL, set to be equal to t.
#' @param centered Boolean, Do the data are centered? (default = FALSE)
#' @param noDiag Boolean, Remove the diagonal of the covariance in the data?
#'  (dafault = TRUE)
#' @param regtype A character string specifying which type of kernel regression
#'  estimator to use. lc specifies a local-constant estimator (Nadaraya-Watson) 
#'  and ll specifies a local-linear estimator. Defaults to lc.
#' 
#' @references Yao, Müller and Wang - Functional Data Analysis for Sparse
#'  Longitudinal Data (2005) - Journal of the American Statistical Association
#' @export
local_covariance <- function(data, t = seq(0, 1, length.out = 101), s = NULL,
                             centered = FALSE, noDiag = TRUE, regtype = 'lc'){
  data_ <- list2cai(data)
  time <- data_$time
  x <- data_$x
  subject <- data_$obs
  
  if (!centered) {
    fit <- np::npreg(x ~ time, ckertype = "epanechnikov", regtype = regtype)
    x <- x - fit$mean
  }
  gg <- NULL
  for(zz in unique(subject)) {
    if(sum(subject == zz) > 1) {
      tt <- time[subject == zz]
      xx <- x[subject == zz]
      g <- expand.grid(t1 = tt, t2 = tt)
      scov <- xx %*% t(xx)
      if(noDiag) 
        scov <- scov + diag(rep(Inf, length(xx)))
      g$z <- matrix(scov, ncol = 1)
      gg <- rbind(gg, g[g$z < Inf, ])
    }
  }

  band <- np::npregbw(xdat = gg[, c(1, 2)], ydat = as.vector(gg[, 3]),
                      ckertype = 'epanechnikov', regtype = regtype)
  
  if(is.null(s))
    s <- t
  new <- expand.grid(t1 = s, t2 = t)
  mod <- np::npreg(band, newdata = new)
  a <- matrix(mod$mean, ncol = length(s), nrow = length(t))
}
