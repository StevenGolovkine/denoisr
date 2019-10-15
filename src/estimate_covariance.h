// estimtate_covariance.h

#ifndef ESTIMATE_COVARIANCE_H
#define ESTIMATE_COVARIANCE_H

#include <RcppArmadillo.h>

arma::mat kernelSmoothingCovariance(
    const List & curves, // Curves list ($x and $t)
    const arma::vec & U, // Estimation points in U
    const arma::vec & V, // Estimation points in V
    const double & b, // Smoothing bandwith for every curve
    const double & h // Global smoothing bandwith
);

#endif