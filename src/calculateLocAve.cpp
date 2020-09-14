
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List calculateLocAve(List nonzeroc, List nonzerog, NumericVector r, List h) {
  Rcpp::List clist(nonzeroc);
  int n = clist.size();
  Rcpp::List glist(nonzerog);
  Rcpp::NumericVector rall(r);
  Rcpp::List hh(h);
  Rcpp::List output(n);
  
  for(int c=0; c<n; c++) {
    SEXP lc = clist[c];
    Rcpp::NumericMatrix cM(lc);
    Rcpp::NumericVector indexc = cM(0,_) - 1;
    Rcpp::NumericVector outc(indexc.size()); //or use cM.ncol()
    for(int g=0; g<indexc.size(); g++){
      SEXP lg = glist[indexc[g]];
      Rcpp::NumericMatrix gM(lg);
      Rcpp::NumericVector indexg = gM(0,_) - 1;
      double ww;
      double outsum = 0;
      double wsum = 0;
      double hg = hh[indexc[g]]; // must declare what type
      for (int k=0; k<indexg.size(); k++) {
        ww = exp(- 0.5 * pow((rall[indexg[k]] - rall[c])/hg, 2) );
        outsum = outsum + (ww * gM(1,k));
        wsum = wsum + ww;
      }
      outc[g] = outsum/wsum;
    }
    output[c] = outc;
  }
  return(output);
}
