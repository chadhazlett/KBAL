#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
// [[Rcpp::depends(RcppParallel)]]
using namespace Rcpp; 
using namespace RcppParallel;


// define both_non_NA(a, b)
inline bool both_non_NA(double a, double b) {
  return (!ISNAN(a) && !ISNAN(b));
}

struct Kernel : public Worker
{
  // source matrix
  const RMatrix<double> X;
  const double b;
  
  // destination matrix
  RMatrix<double> out;
  
  // initialize with source and destination
  Kernel(const Rcpp::NumericMatrix X,
         const double b,
         Rcpp::NumericMatrix out) 
    : X(X), b(b), out(out) {}
  
  // calculate the IBS kernel of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    int p = X.ncol();
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {
        double dist = 0;
        for (int k = 0; k < p; k++) {
          double xi = X(i, k), xj = X(j, k);
          if (both_non_NA(xi, xj)) {
            double diff = xi - xj;
            dist += diff*diff;
          }
        }
        out(i, j) = exp(-dist / b);
        out(j, i) = out(i, j);
      }
      out(i, i) = 1;
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix kernel_parallel(Rcpp::NumericMatrix X,
                                    const double b) {
  
  // allocate the output matrix
  Rcpp::NumericMatrix out(X.nrow(), X.nrow());
  
  // IBSKernel functor (pass input and output matrixes)
  Kernel kernel(X, b, out);
  
  // call parallelFor to do the work
  parallelFor(0, X.nrow(), kernel);
  
  // return the output matrix
  return out;
}

struct Kernel_2 : public Worker
{
  // source matrix
  const RMatrix<double> X;
  const RMatrix<double> Y;
  const double b;
  
  // destination matrix
  RMatrix<double> out;
  
  // initialize with source and destination
  Kernel_2(const Rcpp::NumericMatrix X,
           const Rcpp::NumericMatrix Y,
           const double b,
           Rcpp::NumericMatrix out) 
    : X(X), Y(Y), b(b), out(out) {}
  
  // calculate the IBS kernel of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
      int p = X.ncol();
      int m = Y.nrow();
      for (std::size_t i = begin; i < end; i++) {
          for (std::size_t j = 0; j < m; j++) {
              double dist = 0;
              for (int k = 0; k < p; k++) {
                  double xi = X(i, k), yj = Y(j, k);
                  double diff = xi - yj;
                  dist += diff*diff;
              }
              out(i, j) = exp(-dist / b);
            }
    }
}
};

// [[Rcpp::export]]
Rcpp::NumericMatrix kernel_parallel_2(Rcpp::NumericMatrix X,
                                      Rcpp::NumericMatrix Y,
                                      const double b) {
  
  // allocate the output matrix
  Rcpp::NumericMatrix out(X.nrow(), Y.nrow());
  
  // IBSKernel functor (pass input and output matrixes)
  Kernel_2 kernel_2(X, Y, b, out);
  
  // call parallelFor to do the work
  parallelFor(0, X.nrow(), kernel_2);
  
  // return the output matrix
  return out;
}



struct Kernel_old : public Worker
{
    // source matrix
    const RMatrix<double> X;
    const RMatrix<double> Y;
    const double b;
    
    // destination matrix
    RMatrix<double> out;
    
    // initialize with source and destination
    Kernel_old(const Rcpp::NumericMatrix X,
             const Rcpp::NumericMatrix Y,
             const double b,
             Rcpp::NumericMatrix out) 
        : X(X), Y(Y), b(b), out(out) {}
    
    // calculate the IBS kernel of the range of elements requested
    
    void operator()(std::size_t begin, std::size_t end) {
        int p = X.ncol();
        int m = Y.nrow();
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j < m; j++) {
                double dist = 0;
                for (int k = 0; k < p; k++) {
                    double xi = X(i, k), yj = Y(j, k);
                    if (both_non_NA(xi, yj)) {
                        double diff = xi - yj;
                        dist += diff*diff;
                    }
                }
                out(i, j) = exp(-dist / b);
            }
        }
    }
};


// [[Rcpp::export]]
Rcpp::NumericMatrix kernel_parallel_old(Rcpp::NumericMatrix X,
                                      Rcpp::NumericMatrix Y,
                                      const double b) {
    
    // allocate the output matrix
    Rcpp::NumericMatrix out(X.nrow(), Y.nrow());
    
    // IBSKernel functor (pass input and output matrixes)
    Kernel_old kernel_old(X, Y, b, out);
    
    // call parallelFor to do the work
    parallelFor(0, X.nrow(), kernel_old);
    
    // return the output matrix
    return out;
};

