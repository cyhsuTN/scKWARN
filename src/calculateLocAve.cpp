
#include <Rcpp.h>
using namespace Rcpp;

double weimean(NumericVector x, const double y, const double h, NumericVector z, const double a) {
  NumericVector d1 = (x - y)/h;
  NumericVector d2 = ifelse( abs(x-y) < a * h, exp(- 0.5 * d1 * d1), 0);
  double sumd2 = sum(d2);
  double wm = sum( d2 * z )/sumd2;
  return( wm );
}


// [[Rcpp::export]]
List calculateLocAve(List nonzeroc, List nonzerog, NumericVector r, NumericVector h, const double a) {
  int n = nonzeroc.size();
  Rcpp::List output(n);
  
  for(int c=0; c<n; c++) {
    SEXP lc = nonzeroc[c];
    NumericVector indexc(lc); indexc = indexc - 1;
    NumericVector outc(indexc.size());
    for(int g=0; g<indexc.size(); g++){
      SEXP lg = nonzerog[indexc[g]];
      NumericMatrix gM(lg);
      NumericVector indexg = gM(0,_) - 1;
      double hg = h[indexc[g]];
      outc[g] = weimean(r[indexg], r[c], hg, gM(1,_), a);
    }
    output[c] = outc;
  }
  return(output);
}
