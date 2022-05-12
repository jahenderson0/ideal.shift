############################################
# John Henderson
# - cosineDistF: fast cosine dissimilarity using Rcpp
# - April 18, 2019
############################################

require(Rcpp)
require(RcppArmadillo)
require(inline)

############################################
#FAST COSINE DISSIMILARITY FUNCTION -- to use in preliminary dimension reduction
############################################

fastdistC='
  Rcpp::NumericMatrix dtm(wfm);
  int N = dtm.nrow();
  int M = dtm.ncol();
  Rcpp::NumericMatrix dmat(N,N);
  Rcpp::NumericVector rlen(N);
  arma::mat ip;
  arma::mat a;
  arma::mat b;
  for (int i=0; i < N; i++){
    rlen(i)=sqrt(sum(pow(dtm(i,_),2)));
  }
  for (int i=0; i < (N-1); i++){
    a = dtm(i,_);
    for (int j=(i+1); j < N; j++){
      b = dtm(j,_);
      ip=a.t()*b;
      dmat(i,j) = 1-ip(0)/(rlen(i) * rlen(j));
      dmat(j,i) = dmat(i,j);
    }
  }

  return Rcpp::List::create(Rcpp::Named("dmat") = dmat);
'

cosineDistF = cxxfunction(signature(wfm = "numeric"), body=fastdistC, plugin = "RcppArmadillo")

#END cosineDistF
