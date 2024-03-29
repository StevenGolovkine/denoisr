% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_H0.R
\name{estimate_H0_deriv}
\alias{estimate_H0_deriv}
\title{Perform an estimation of \eqn{H_0} when the curves are derivables}
\usage{
estimate_H0_deriv(data, t0 = 0, eps = 0.01, k0 = 2, sigma = NULL)
}
\arguments{
\item{data}{A list, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} The sampling points
\item \strong{$x} The observed points.
}}

\item{t0}{Numeric, the sampling points at which we estimate \eqn{H_0}. We will
consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for the estimation of
\eqn{H_0} when \eqn{\sigma} is unknown.}

\item{eps}{Numeric, precision parameter. It is used to control how much larger
than 1, we have to be in order to consider to have a regularity larger than 1
(default to 0.01).}

\item{k0}{Numeric, the number of neighbors of \eqn{t_0} to consider.}

\item{sigma}{Numeric, true value of sigma. Can be NULL.}
}
\value{
Numeric, an estimation of \eqn{H_0}.
}
\description{
This function performs an estimation of \eqn{H_0} used for the estimation of
the bandwidth for a univariate kernel regression estimator defined over
continuous domains data in the case the curves are derivables.
}
\examples{
X <- generate_integrate_fractional_brownian(N = 1000, M = 300,
                                            H = 0.5, sigma = 0.01)
estimate_H0_deriv(X, t0 = 0.5, eps = 0.01, k0 = 6)
}
\references{
Golovkine S., Klutchnikoff N., Patilea V. (2020) - Learning the
smoothness of noisy curves with applications to online curves denoising.
}
\seealso{
Other estimate \eqn{H_0}: 
\code{\link{estimate_H0_deriv_list}()},
\code{\link{estimate_H0_list}()},
\code{\link{estimate_H0}()}
}
\concept{estimate \eqn{H_0}}
