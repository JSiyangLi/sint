#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void iterateCombinations(IntegerVector max_vals, Function process_combo) {
  int n = max_vals.size();
  IntegerVector current(n, 0);  // Start at [0, 0, ..., 0]
  
  while (true) {
    // Call R function with current combo
    process_combo(current);
    
    // Move to next combination
    bool exhausted = true;
    for (int i = n - 1; i >= 0; --i) {
      if (current[i] < max_vals[i]) {
        current[i]++;
        std::fill(current.begin() + i + 1, current.end(), 0);
        exhausted = false;
        break;
      }
    }
    
    if (exhausted) break;
  }
}