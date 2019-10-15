// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_covariance.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "estimate_curve.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat kernelSmoothingCovariance(
    const List & curves, // Curves list ($x and $t)
    const arma::vec & U, // Estimation points in U
    const arma::vec & V, // Estimation points in V
    const double & b, // Smoothing bandwith for every curve
    const double & h // Global smoothing bandwith
  ){

  // Get parameters
  arma::uword N = curves.length(); // Number of curves

  // Create a vector containing all the sampling points (unique)
  List mycurve = curves[0];
  arma::vec T = as<arma::vec>(mycurve["t"]);

  Rcout << "=== Creation of the sampling vector ===" << std::endl;
  for(arma::uword n=1; n<N; n++){
    Rcout << (n+1)*100/N << "%  \r";

    mycurve = curves[n]; // Get the nth curve
    T = arma::join_cols(T, as<arma::vec>(mycurve["t"]));
  }
  T = arma::unique(T); // unique() returns sorted elements
  Rcout << "\n=== End ===" << std::endl;
  
  // Smooth each curve 
  arma::mat Y_hat(N, T.n_elem); Y_hat.fill(0);

  Rcout << "=== Smoothing of the curves ===" << std::endl;
  for(arma::uword n=0; n<N; n++){
    Rcout << (n+1)*100/N << "% \r";

    mycurve = curves[n];
    arma::vec x = mycurve["x"];
    arma::vec t = mycurve["t"];

    Y_hat.row(n) = kernelSmoothingCurve(T, t, x, b).t();
  }
  Rcout << "\n=== End ===" << std::endl;

  // Compute the leave-one-out means
  arma::mat Y_bar_LOO(N, T.n_elem); Y_bar_LOO.fill(0);
  arma::mat Y_bar = arma::mean(Y_hat, 0); // (1 x T.n_elem) matrix 

  Rcout << "=== Leave-one-out means ===" << std::endl;
  for(arma::uword n=0; n<N; n++){
    Rcout << (n+1)*100/N << "% \r";
    Y_bar_LOO.row(n) = N * Y_bar / (N - 1) - Y_hat.row(n) / (N - 1);
  }
  Rcout <<"\n=== End ===" << std::endl;
  
  // Get the means at the right index
  List Y_bbar(N);

  Rcout << "=== Indexing the means ===" << std::endl;
  for(arma::uword n=0; n<N; n++){
    Rcout << (n+1)*100/N << "% \r";

    mycurve = curves[n];
    arma::vec t = mycurve["t"];
    arma::vec v(t.n_elem); v.fill(0);

    for(arma::uword i=0; i<t.n_elem; i++){
      arma::uvec idx = find(T == t(i));
      v(i) = Y_bar_LOO(n, idx(0));
    }
    Y_bbar[n] = v.col(0);
  }
  Rcout << "\n=== End ===" << std::endl;

  // Estimate covariance
  arma::mat phi_hat(U.n_elem, V.n_elem); phi_hat.fill(0);
  double kernel = 0;
  Rcout << "=== Covariance estimation ===" << std::endl;
  for(arma::uword iu=0; iu<U.n_elem; iu++){
    Rcout << (iu+1)*100/U.n_elem << "% \r";
    
    for(arma::uword iv=0; iv<V.n_elem; iv++){
      
      double inter2 = 0;
      for(arma::uword n=0; n<N; n++){

        mycurve = curves[n];
        arma::vec X_n = mycurve["x"];
        arma::vec T_n = mycurve["t"];

        arma::uword M_n = X_n.n_elem;

        arma::vec mean_curve = Y_bbar[n];

        double inter1 = 0;
        for(arma::uword k=0; k<M_n; k++){
          double K_k_u = 0;
          if(std::abs(T_n(k) - U(iu)) < h){
            K_k_u = 3 * (1 - std::pow((T_n(k) - U(iu))/h, 2)) / (4 * h);
          }

          arma::vec Y_hat_kl(2); Y_hat_kl.fill(0);
          for(int l=0; l<M_n; l++){
            if(l != k){

              // arma::vec XX(M_n - 2); XX.fill(0);
              // arma::vec TT(M_n - 2); TT.fill(0);
              // arma::uvec idx(M_n - 2); idx.fill(0);
              // arma::uword I = 0;
              // for(arma::uword i=0; i<M_n; i++){
              //   if(i != k && i != l){
              //     idx(I++) = i;
              //   }
              // }
              // XX = X_n.elem(idx);
              // TT = T_n.elem(idx);

              Y_hat_kl = kernelSmoothingCurve({T_n(k), T_n(l)}, T_n, X_n, b);

              double K_l_v = 0;
              if(std::abs(T_n(l) - V(iv)) < h){
                K_l_v = 3 * (1 - std::pow((T_n(l) - V(iv))/h, 2)) / (4 * h);
              }
              
              inter1 += (Y_hat_kl[0] - mean_curve[k]) * (Y_hat_kl[1] - mean_curve[l]) * K_k_u * K_l_v;
              kernel += K_k_u * K_l_v;
            }
          }

        }
        inter2 += (inter1 / (M_n * (M_n - 1)));
      }
      phi_hat(iu, iv) = (inter2 / N);
    }
  }
  Rcout << "\n=== End ===" << std::endl;

  return(phi_hat);
}