// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// mfBm.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>


using namespace Rcpp;

// [[Rcpp::export]]
arma::vec mfBm(
    arma::vec & H, 
    arma::vec & T,
    int N,
    int K) {
  
  arma::vec res(K);
  
  double tmp1, tmp2, tmp3, tmp;
  double R(0);
  std::srand(time(NULL));
  
  /* Loop over the terms of the sum (de 0 to N-1) */
  for(arma::uword n=0; n < N; n++){
    
    double theta = ((double) std::rand() / (RAND_MAX));
    double U(1);
    if ((rand()/(float)RAND_MAX) < 0.5){
      U = -1;
    }
    R = R - log(rand()/(float)RAND_MAX);
    
    /* Loop over each discretization points (0 à K-1) */
    for(arma::uword k=0; k < K; k++){
      tmp1 = cos(2 * M_PI * theta) * (cos(T[k] * R * U) - 1);
      tmp2 = sin(2 * M_PI * theta)*sin(-T[k] * R * U);
      tmp3 = pow(R,0.5 + H[k]);
      tmp = (tmp1 - tmp2) / tmp3;
      res[k] = res[k] + tmp;
      if (n == N - 1){
        res[k] = res[k] * 2;
      }
    } 
  }
    
  return(res);
}

/*** R
plot(mfBm(rep(0.5, 1000), seq(from=0,to=1,length=1000), 1000, 50))

*/
