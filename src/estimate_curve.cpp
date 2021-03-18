// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_curve.cpp
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec epaKernelSmoothingCurve(
    const arma::vec & U, // Estimation points in U 
    const arma::vec & T, // Sampling points
    const arma::vec & Y, // Curves points
    const arma::vec & b // Smoothing bandwiths
  ){

  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double cpt = 0;
    // Loop over the known points
    for (arma::uword k=0; k<M; k++){
      if (std::abs(T(k) - U(i)) <= b(i)){
        //if (std::abs(T(k) - U(i)) != 0){
        Y_hat(i) += (Y(k) * (1 - std::pow((T(k) - U(i))/b(i), 2)));
        cpt += (1 - std::pow((T(k) - U(i))/b(i), 2));
        //}
      }
    }
    if(cpt == 0){
      Y_hat(i) = R_NaN;
    } else{
      Y_hat(i) /= cpt;
    }
  }
  
  return Y_hat;
}

// [[Rcpp::export]]
arma::vec uniKernelSmoothingCurve(
    const arma::vec & U, // Estimation points in U 
    const arma::vec & T, // Sampling points
    const arma::vec & Y, // Curves points
    const arma::vec & b // Smoothing bandwiths
){
  
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double cpt = 0;
    // Loop over the known points
    for (arma::uword k=0; k<M; k++){
      
      if (std::abs(T(k) - U(i)) <= b(i)){
        Y_hat(i) += Y(k);
        cpt += 1;
      }
    }
    if(cpt == 0){
      Y_hat(i) = R_NaN;
    } else{
      Y_hat(i) /= cpt;
    }
  }
  
  return Y_hat;
}


// [[Rcpp::export]]
arma::vec betaKernelSmoothingCurve(
  const arma::vec & U, // Estimation points 
  const arma::vec & T, // Sampling points
  const arma::vec & Y, // Curves points
  const arma::vec & b // Smoothing bandwiths
){
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Create the S vector
  arma::vec S(M);
  S.fill(0); S(M-1) = 1;
  for(arma::uword m=1; m<(M-1); m++){
    S(m) = (T(m) + T(m+1)) / 2;
  }
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double B = R::beta(U(i)/b(i) + 1, (1 - U(i))/b(i) + 1);
    // Loop over the known points
    for(arma::uword k=1; k<M; k++){
      
      arma::vec S_k = arma::linspace<arma::vec>(S(k-1), S(k), 100);
      arma::vec K = arma::pow(S_k, U(i)/b(i)) % arma::pow((1 - S_k), (1 - U(i))/b(i));
      arma::mat intK = arma::trapz(S_k, K);
      
      Y_hat(i) += (Y(k) * intK(0, 0)) / B;
    }
  }
  
  return Y_hat;
}
  
// [[Rcpp::export]]
arma::vec modifiedBetaKernelSmoothingCurve(
  const arma::vec & U, // Estimation points
  const arma::vec & T, // Sampling points
  const arma::vec & Y, // Curves points
  const arma::vec & b // Smoothing bandwiths
){
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points 
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Create the S vector
  arma::vec S(M);
  S.fill(0); S(M-1) = 1;
  for(arma::uword m=1; m<(M-1); m++){
    S(m) = (T(m) + T(m+1)) / 2;
  }
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    
    if(U(i) < 2*b(i)){
      double rho = 2 * std::pow(b(i), 2) + 2.5 - std::pow(4 * std::pow(b(i), 4) + 6 * std::pow(b(i), 2) + 2.25 - std::pow(U(i), 2) - U(i) / b(i), 0.5);
      double B = R::beta(rho, (1 - U(i))/b(i));
      // Loop over the known points
      for(arma::uword k=1; k<M; k++){
        
        arma::vec S_k = arma::linspace<arma::vec>(S(k-1), S(k), 10);
        arma::vec K = arma::pow(S_k, rho - 1) % arma::pow((1 - S_k), (1 - U(i))/b(i) - 1);
        arma::mat intK = arma::trapz(S_k, K);
        
        Y_hat(i) += (Y(k) * intK(0, 0)) / B;
      }
    } else if(U(i) > 1 - 2*b(i)){
      double rho = 2 * std::pow(b(i), 2) + 2.5 - std::pow(4 * std::pow(b(i), 4) + 6 * std::pow(b(i), 2) + 2.25 - std::pow((1 - U(i)), 2) - (1 - U(i)) / b(i), 0.5);
      double B = R::beta(U(i)/b(i), rho);
      // Loop over the known points
      for(arma::uword k=1; k<M; k++){
        
        arma::vec S_k = arma::linspace<arma::vec>(S(k-1), S(k), 10);
        arma::vec K = arma::pow(S_k, U(i)/b(i) - 1) % arma::pow((1 - S_k), rho - 1);
        arma::mat intK = arma::trapz(S_k, K);
        
        Y_hat(i) += (Y(k) * intK(0, 0)) / B;
      }
    } else{ // 2*b < U(i) < 1 - 2*b
      double B = R::beta(U(i)/b(i), (1 - U(i))/b(i));
      // Loop over the known points
      for(arma::uword k=1; k<M; k++){
        
        arma::vec S_k = arma::linspace<arma::vec>(S(k-1), S(k), 10);
        arma::vec K = arma::pow(S_k, U(i)/b(i) - 1) % arma::pow((1 - S_k), (1 - U(i))/b(i) - 1);
        arma::mat intK = arma::trapz(S_k, K);
        
        Y_hat(i) += (Y(k) * intK(0, 0)) / B;
      }
    }

  }
  
  return Y_hat;
}
