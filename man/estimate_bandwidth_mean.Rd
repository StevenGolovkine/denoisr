% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_b.R
\name{estimate_bandwidth_mean}
\alias{estimate_bandwidth_mean}
\title{Perform an estimation of the bandwidth for the estimation the mean}
\usage{
estimate_bandwidth_mean(
  data,
  t0_list = 0.5,
  k0_list = 2,
  grid = NULL,
  nb_obs_minimal = 2,
  type_k = 2
)
}
\arguments{
\item{data}{A list, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} The sampling points
\item \strong{$x} The observed points.
}}

\item{t0_list}{A vector of numerics, the sampling points at which we estimate
\eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for
the estimation of \eqn{L_0} when \eqn{\sigma} is unknown.}

\item{k0_list}{A vector of numerics, the number of neighbors of \eqn{t_0} to
consider. Should be set as \deqn{k0 = (M / log(M) + 7) / 8}. We can set a
different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just
put a unique numeric.}

\item{type_k}{Integer, used kernel!. If 1, uniform kernel, If 2, Epanechnikov
kernel. If 3, biweight kernel.}
}
\value{
A list, with elements:
\itemize{
\item \strong{sigma} An estimation of the standard deviation of the noise
\item \strong{H0} An estimation of \eqn{H_0}
\item \strong{L0} An estimation of \eqn{L_0}
\item \strong{b} An estimation of the bandwidth
}
}
\description{
This function performs an estimation of the bandwidth to be used in the
Nadaraya-Watson estimator. The bandwidth is estimated using the method from
\cite{add ref}.
}
