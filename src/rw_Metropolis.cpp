#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' @name rw_Metropolis
//' @title A Metropolis sampler using Rcpp
//' @description A Metropolis sampler using Rcpp
//' @param N the number of samples
//' @param x0 the number of between-sample random number
//' @param sigma variance
//' @return a List with a random sample of size \code{n}
//' @export
//[[Rcpp::export]]
double dLaplace(double x) {
  return 0.5*exp(-abs(x));
}

//[[Rcpp::export]]
List rw_Metropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0; 
  int k=0;
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (dLaplace(y[0]) / dLaplace(x[i-1]))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1]; 
      k++;
    }
  }
  return List::create(
    _["x"] = x,
    _["k"] = k
  );
} 
