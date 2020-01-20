################################################################################
#            Functions for risk estimation using regularity                    #
################################################################################


#' Perform an estimation of the risk on a set of curves at a particular point.
#' 
#' This function performs the estimation of the risk on a set of curves at a
#' particular sampling points \eqn{t_0}. Both the real and estimated curves have
#' to be sampled on the same grid. 
#' 
#' Actually, two risks are computed. They are defined as:
#' \deqn{MeanRSE(t_0) = \frac{1}{N}\sum_{n = 1}^{N}(X_n(t_0) - \hat{X}_n(t_0))^2}
#' and
#' \deqn{MeanRSE(t_0) = \max_{1 \leq n \leq N} (X_n(t_0) - \hat{X}_n(t_0))^2}
#'
#' @param curves A list, where each element represents a real curve. Each curve 
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param curves_estim A list, where each element represents an estimated curve. 
#'  Each curve have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  } 
#' @param t0_list A vector of numerics, sampling points where the risk will be 
#'  computed. Can have a single value.
#' 
#' @return A list, with the mean and max residual squared error in \eqn{t_0}.
#' @export
#' @examples 
#' data("fractional_brownian")
#' curves_smoothed <- smooth_curves(fractional_brownian)$smooth
#' estimate_risk(fractional_brownian, curves_smoothed, t0_list = 0.5)
estimate_risk <- function(curves, curves_estim, t0_list = 0.5) {
  
  risk_df <- dplyr::tibble(t0 = numeric(), 
                           MeanRSE = numeric(), 
                           MaxRSE = numeric())
  
  for(t0 in t0_list) {
    risk <- estimateRisk(curves, curves_estim, t0)
    risk_df <- risk_df %>% dplyr::add_row(t0 = t0, MeanRSE = risk[1], MaxRSE = risk[2])
  }

  risk_df
}

#' Perform an estimation of the risk on a set of curves along the sampling points.
#' 
#' This function performs the estimation of the risk on a set of curves along 
#' the sampling points. Both the real and estimated curves have to be sampled on
#' the same grid.
#' 
#' Actually, two risks are computed. They are defined as:
#' \deqn{MeanIntRSE = \frac{1}{N}\sum_{n = 1}^{N}\int(X_n(t) - \hat{X}_n(t))^2dt}
#' and
#' \deqn{MeanIntRSE = \max_{1 \leq n \leq N} \int(X_n(t) - \hat{X}_n(t))^2dt}
#'
#' @param curves A list, where each element represents a real curve. Each curve 
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param curves_estim A list, where each element represents an estimated curve. 
#'  Each curve have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  }
#'
#' @return A list, with the mean and max integrated residual squared error 
#'  in \eqn{t_0}.
#' @export
#' @examples 
#' data("fractional_brownian")
#' curves_smoothed <- smooth_curves(fractional_brownian)$smooth
#' estimate_risks(fractional_brownian, curves_smoothed)
estimate_risks <- function(curves, curves_estim) {
  risk <- estimateRiskCurves(curves, curves_estim)

  c("MeanIntRSE" = risk[1], "MaxIntRSE" = risk[2])
}


#' Perform an estimation of the risk on one curve along the sampling points.
#' 
#' This function performs the estimation of the risk on one curve along 
#' the sampling points. Both the real and estimated curve have to be sampled on
#' the same grid.
#' 
#' Actually, one risk is computed. They are defined as:
#' \deqn{IntRSE = \int(X_n(t) - \hat{X}_n(t))^2dt}
#' 
#' @param curves A list, where each element represents a real curve. Each curve 
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param curves_estim A list, where each element represents an estimated curve. 
#'  Each curve have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  }
#'
#' @return Numeric, the integrated mean squarred eror
#' @export
#' @examples 
#' data("fractional_brownian")
#' curves_smoothed <- smooth_curves(fractional_brownian)$smooth
#' estimate_risks(fractional_brownian[[1]], curves_smoothed[[1]])
estimate_int_risk <- function(curves, curves_estim) {
  estimateRiskCurve(curves, curves_estim)
}

