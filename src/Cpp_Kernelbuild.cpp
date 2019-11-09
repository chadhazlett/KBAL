#include <Rcpp.h>
using namespace Rcpp;

/*
double euc_dist(const arma::rowvec& x1, const arma::rowvec& x2) {
  double out = 0.0;
  unsigned n = x1.n_elem;

  for (unsigned i = 0; i < n; ++i) {
    out += pow(x1(i) - x2(i), 2);
  }
  return sqrt(out);
}

double kern_gauss_1d(const arma::rowvec& x1, const arma::rowvec& x2, const double& b)
{
  return exp(-pow(euc_dist(x1, x2), 2) / (b));
}
*/

//[[Rcpp::export]]
NumericMatrix new_gauss_kern(const NumericMatrix newx, const NumericMatrix oldx, const double b) {
  unsigned n1 = newx.nrow();
  unsigned n2 = oldx.nrow();
  NumericMatrix out(n1, n2);
  double euc_distsq = 0;
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      for(int k = 0; k < newx.ncol(); ++k) {
          //Rcout << newx(i,k) - oldx(j,k) << std::endl; uncomment for printout debugging
          euc_distsq += pow( (newx(i,k) - oldx(j,k)), 2.0);
          //ALT:no separate k loop ad use vectors subtracted? is that faster?
      }
      out(i,j) = exp(-euc_distsq/b);
      euc_distsq = 0;
    }

  }
  return out;
}

/* to test internally when compile source code uncoment all of this and 
  adjust the *s accordingly so there are no spaces anymore
# / * ** R
# data(lalonde)
# lalonde$nodegr=as.numeric(lalonde$educ<=11)
# xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
# allx = lalonde[,xvars]
# #allx_mini <- allx[1:10,]
# treatment=lalonde$nsw[]
# observed = 1-treatment
# target = treatment
# useasbases=as.numeric(observed==1)
# b=2*ncol(allx)
# 
# bases = allx[useasbases==1, ]
# X = bases
# newData = allx
# 
# Xmeans <- colMeans(X)
# Xsd <- apply(X, 2, sd)
# X <- scale(X, center = Xmeans, scale = Xsd)
# newData <- scale(newData, center = Xmeans, scale = Xsd)
# 
# newK <- new_gauss_kern(newData, X, b)
# * /
*/
