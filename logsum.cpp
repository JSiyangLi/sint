//#include <Rcpp.h>
//using namespace Rcpp;
#include <iostream>
using namespace std;

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
double logplusvec(vector<double> & x){
  
  int n = x.size();
  
  double r = -numeric_limits<double>::max();
  
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

// python version
extern "C" {
    double py_logplus(double x, double y){
        return logplus(x, y);
    }

    double py_logplusvec(vector<double> & x){
        return logplusvec(x);
    }

    double py_logminus(double x, double y){
        return logminus(x, y);
    }
}
