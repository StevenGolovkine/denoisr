################################################################################
#                  Functions for covariance estimation with RKHS               #
################################################################################

#' Estimation of the mean with spline smoothing
#' 
#' This function performs the estimation of the mean function based on the paper
#' by Cai and Yuan who show the optimal estimation of the mean function under
#' both common and independent designs.
#' 
#' @param data A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param time A vector containing the times where the mean is computed.
#' 
#' @return A vector representing the mean curve.
#'  
#' @references Cai and Yuan, Optimal Estimation of the mean function based on
#'  discretely sampled functional data: phase transition, 2011, 
#'  The Annals of Statistics
#' @export
ssf_mean <- function(data, time){
  data_ <- list2cai(data)
  mod <- stats::smooth.spline(data_$time, data_$x)
  stats::predict(mod, time)
}

#' Estimation of the covariance with spline smoothing
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
#' @param nbasis Integer, number of basis to use for the splines.
#' @param centered Boolean, Do the data are centered? (default = FALSE)
#' @param noDiag Boolean, Remove the diagonal of the covariance in the data?
#'  (dafault = TRUE)
#' 
#' @references Cai and Yuan, Nonparametric Covariance Function Estimation for
#' Functional and Longitudinal Data, 2010.
#' @export
ssf_covariance <- function(data, t = seq(0, 1, length.out = 101), s = NULL,
                           nbasis = 5, centered = FALSE, noDiag = TRUE){
  predict_ssanova <- utils::getFromNamespace("predict.ssanova", "gss")
  
  data_ <- list2cai(data)
  time <- data_$time
  x <- data_$x
  subject <- data_$obs
  
  if (!centered) {
    fit <- stats::smooth.spline(time, x)
    x <- x - stats::fitted(fit)
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
  
  gg <- unique(gg)
  nobs <- nrow(gg)
  tt <- min(time) + (max(time) - min(time)) * (1:nbasis)/(nbasis + 1)
  g <- expand.grid(t1 = tt, t2 = tt)
  g$z <- 0
  gg <- rbind(g, gg)
  
  fit <- gss::ssanova(z ~ t1 * t2, data = gg, id.basis = 1:(nbasis * nbasis))
  
  if(is.null(s))
    s <- t
  new <- expand.grid(t1 = s, t2 = t)
  estim <- predict_ssanova(fit, newdata = new)
  matrix(estim, ncol = length(s), nrow = length(t))
}
