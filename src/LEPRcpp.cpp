#define ARMA_DONT_USE_WRAPPER
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "Rcpp_aux.hpp"
#include "aux.hpp"
#include "method.hpp"
#include "readPlink.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP LEP(SEXP Xs, arma::vec& y, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12, bool verbose = false, bool memory_copy = true){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int method = LEP_METHOD;
  int N, P;
  int type;
  void* ptr;
  List ret;
  Col<uword> xindex;
  preprocess_summary(Xs, SS, lbPval, type, N, P, summary, ptr, xindex, verbose);
  FitObj* fitobj = NULL;
  if(Rf_isInteger(Xs)){
    int* data = INTEGER(Xs);
    fitobj = fit(data, y, P, &summary, &xindex, lp_opt, type, verbose, method);
  }else if(Rf_isReal(Xs)){
    double* data = NULL;
    if(!memory_copy){
      data = REAL(Xs);
    }else{
      data = new double[(unsigned long long)N * (unsigned long long)P];
      memcpy(data, REAL(Xs), sizeof(double) * (unsigned long long)N * (unsigned long long)P);
    }
    fitobj = fit(data, y, P, &summary, &xindex, lp_opt, type, verbose, method);
    if(memory_copy){
      delete[] data;
    }
  }else{
    cout << "Wrong type of input, the input matrix should be double precision matrix or integer matrix!" << endl;
    return ret;
  }

  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  ret = wrap_fit(fitobj);
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP LEPCV(SEXP Xs, arma::vec& y, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12, Rcpp::String measure = "mse", bool verbose = false){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int method = LEP_METHOD;
  int N, P, type;
  void* ptr;
  List ret;
  Col<uword> xindex;
  preprocess_summary(Xs, SS, lbPval, type, N, P, summary, ptr, xindex, verbose);
  double predict = 0;
  if(Rf_isInteger(Xs)){
    int* data = INTEGER(Xs);//static_cast<Mat<double> *>(ptr) ->memptr();
    predict = CV_Proc(data , y, P, &summary, &xindex, lp_opt, type, measure, verbose, method);
  }else if(Rf_isReal(Xs)){
    double* data = REAL(Xs);
    predict = CV_Proc(data , y, P, &summary, &xindex, lp_opt, type, measure, verbose, method);
  }else{
    cout << "Wrong type of input, the input matrix should be double precision matrix or integer matrix!" << endl;
    return ret;
  }

  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  ret[measure] = wrap(predict);
  return ret;
}

// [[Rcpp::export]]
RcppExport SEXP LEP_Plink(Rcpp::String genoplinkfile, SEXP  SS = R_NilValue, SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12, bool verbose = false){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int method = LEP_METHOD;
  int N, P;
  int type = type_char;
  List ret;
  vec y;
  Chroms chroms;
  chroms = getChromsY(genoplinkfile, y, N, P);
  long long size = (long long)N * (long long)P;
  char* X = new char[ size];
  clock_t t1 = clock();
  readPlink(genoplinkfile,N, P, X);
  if(verbose){
    cout <<"Time for reading the data is " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
  }
  Col<uword> xindex;
  summary = combine_summary_X(chroms.snps, SS, lbPval, verbose, P, xindex);
  // cout << "chroms.snps = " << chroms.snps.size() << " " << " xindex = " << xindex.size() << endl;
  if( verbose ){
    cout << "dimension of summaries : " << summary.n_rows << " " << summary.n_cols << endl;
    cout << "intersects of genotype data and summary data is " << xindex.size() << endl;
  }
  FitObj* fitobj = fit(X, y, P, &summary, &xindex,lp_opt, type, verbose, method);
  vector<string> snps_names;
  if(xindex.size() > 0){
    for(int i = 0; i < xindex.size(); i++){
      snps_names.push_back(chroms.snps.at(xindex.at(i)).name);
    }
  }else{
    for(int i = 0; i < chroms.snps.size(); i++){
      snps_names.push_back(chroms.snps.at(i).name);
    }
  }


  ret = wrap_fit(fitobj);
  ret["snpnames"] = snps_names;
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  delete[] X;
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP LEPCV_Plink(Rcpp::String genoplinkfile, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12,  Rcpp::String measure = "mse", bool verbose = false){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  int method = LEP_METHOD;
  mat summary;
  int N, P;
  int type = type_char;
  List ret;
  vec y;
  Chroms chroms;
  chroms = getChromsY(genoplinkfile, y, N, P);
  long long size = (long long)N * (long long)P;
  char* X = new char[ size];
  clock_t t1 = clock();
  readPlink(genoplinkfile,N, P, X);
  cout <<"Time for reading the data is " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
  Col<uword> xindex;
  summary = combine_summary_X(chroms.snps, SS, lbPval, verbose, P, xindex);
  double predict =CV_Proc(X, y, P, &summary, &xindex, lp_opt, type, measure, verbose, method);
  ret[measure] = wrap(predict);
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  delete[] X;
  return ret;
}



// [[Rcpp::export]]
RcppExport SEXP Predict(SEXP fit_,  arma::mat& X){
  vec ypred(X.n_rows);
  FitObj* fit = NULL;
  if(!Rf_isNull(fit_)){
    Rcpp::List fitList(fit_);
    fit = new FitObj(fitList["gammas"], fitList["mu"], fitList["cov"]);
    ypred = fit -> predict(&X);
  }else{
    cout << "Invalid input of LEP fit!" << endl;
    return Rcpp::wrap(ypred);
  }
  return Rcpp::wrap(ypred);
}
