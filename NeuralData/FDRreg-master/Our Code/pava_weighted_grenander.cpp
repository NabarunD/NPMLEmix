#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector droplast(NumericVector v){
  int l = v.size();
  NumericVector vnew = v[Range(0,l-2)];
  return vnew;
}

// [[Rcpp::export]]
NumericVector pwg(NumericVector w, NumericVector d){

  int n = d.size();
  NumericVector fxx = w/d;
  fxx = fxx/sum(fxx * d);
  
  int num_blocks = n;
  NumericVector block_counts(n);
  for (int i = 0; i < n; i++)
    block_counts[i] = 1;
  NumericVector block_cumsum(n);
  for (int i = 0; i < n; i++){
    block_cumsum[i] = i + 1;
  }

  bool flag = TRUE;
  while(flag){
    // find left-most block where f is locally increasing
    int badblock = 0;
    bool flag_findblock = TRUE;
    while(flag_findblock){
      if (fxx[badblock] < fxx[badblock + 1]) {
        flag_findblock = FALSE;
      } else {
        badblock++;
        if (badblock == num_blocks - 1){
          flag_findblock = FALSE;
        }
      }
    }
    // finished finding badblock
    
    // writing some output functions
    // Rcout << "block count is" << std::endl << block_counts << std::endl;
    // Rcout << "block cumsum is" << std::endl << block_cumsum << std::endl;
    // Rcout << "fxx is" << std::endl << fxx << std::endl;
    // Rcout << "badblock is" << std::endl << badblock << std::endl;
    // finished outputting

    if(badblock < num_blocks - 1){
      num_blocks--;
      // merge the two blocks situated at badblock and badblock + 1
      // updating block counts and cumsums
      block_counts[badblock] += block_counts[badblock+1];
      block_cumsum[badblock] = block_cumsum[badblock+1];
      for (int i = badblock + 1; i < num_blocks; i++){
        block_counts[i] = block_counts[i+1];
        block_cumsum[i] = block_cumsum[i+1];
      }

      // compting the new value of f on this block
      int blower = block_cumsum[badblock] - block_counts[badblock];
      int bupper = block_cumsum[badblock] - 1;
      double bsw = sum(w[Range(blower, bupper)]);
      double bsd = sum(d[Range(blower, bupper)]);
      double newfb = bsw/bsd;

      fxx[badblock] = newfb;
      for (int i = badblock + 1; i < num_blocks; i++){
        fxx[i] = fxx[i+1];
      }

    } else {
      flag = FALSE;
    }
  }

  NumericVector fxxout(n);
  for (int i = 0; i < num_blocks; i++){
    int blower = block_cumsum[i] - block_counts[i];
    int bupper = block_cumsum[i] - 1;
    for (int j = blower; j <= bupper; j++){
      fxxout[j] = fxx[i];
    }  
  }

  return fxxout;
}
