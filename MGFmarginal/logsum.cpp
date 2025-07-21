#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double logplus(double x, double y){
  
  if(x>y){
    
    return x + log(1.0+exp(y-x));
    
  }else{
    
    return y + log(1.0+exp(x-y));
    
  }
  
}
// [[Rcpp::export]]
double logplusvec( NumericVector & x){
  
  int n = x.size();
  
  double r = -DBL_MAX;
  
  for(int i = 0; i < n; i++){
    
    r = logplus(r, x[i]);
    
  }
  
  return r;
  
}
// [[Rcpp::export]]
double logminus(double x, double y){
  
  if(x>=y){
    
    return x + log(1.0-exp(y-x));
    
  }else{
    
    return NAN;
    
  }
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

