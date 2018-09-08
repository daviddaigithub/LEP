#ifndef Rcpp_aux_hpp
#define Rcpp_aux_hpp
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "aux.hpp"
#include "readPlink.hpp"
using namespace arma;
using namespace Rcpp;
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/algorithm/string.hpp>


int snps_overlap(vector<SNP>& x_snps, CharacterVector& ss_snps, Col<uword>& xindex, Col<uword>& sindex);


Options* read_opts(SEXP opts){
  Options* lp_opt = NULL;
  if(!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    int n_fold = 5;
    int max_iter = 300;
    int display_gap= 60;
    if(opt.containsElementNamed("n_fold")){
      n_fold = opt["n_fold"];
    }
    if(opt.containsElementNamed("max_iter")){
      max_iter = opt["max_iter"];
    }

    if(opt.containsElementNamed("dis_gap")){
      display_gap = opt["dis_gap"];
    }
    lp_opt = new Options(max_iter, display_gap, n_fold);
  }
  return lp_opt;

}

void getXNames(SEXP Xs, vector<string>& xnames){
  if(!Rf_isNull(Xs) && !Rf_isMatrix(Xs)){
    Function asMatrix("as.matrix");
    Xs = asMatrix(Xs);
  }
  CharacterVector geno_snps = colnames(Xs);
  for(int i = 0; i < geno_snps.length(); i++){
    xnames.push_back(std::string(geno_snps[i]));
  }

}


mat getSummary(SEXP SS, double lbPVal, vector<string>& snp_names){
  CharacterVector snps_from_ss;
  if(!Rf_isNull(SS) && !Rf_isMatrix(SS)){
    Function asMatrix("as.matrix");
    SS = asMatrix(SS);
  }
  if(!Rf_isNull(SS) ){
    if(!Rf_isNull(rownames(SS))){
      snps_from_ss = rownames(SS);
      for(int i = 0; i < snps_from_ss.length(); i++){
        snp_names.push_back(std::string(snps_from_ss[i]));
      }
    }
  }
  mat lp_summaries(0,0);
  if(!Rf_isNull(SS)){
    NumericMatrix SS_Matrix(SS);
    mat lp_summaries(SS_Matrix.begin(), SS_Matrix.nrow(), SS_Matrix.ncol(), false);
    if(lp_summaries.n_cols > 0 && lp_summaries.n_rows > 0){
      uvec q = find(lp_summaries < lbPVal);
      lp_summaries.elem(q).fill(lbPVal);
    }
    if(lp_summaries.n_rows > lp_summaries.n_cols){
      lp_summaries = lp_summaries.t();
    }
    return lp_summaries;
  }
  return lp_summaries;


}


mat combine_summary_X(vector<SNP>& X_snps, SEXP SS, double lbPval, bool verbose, int P, Col<uword>& xindex)
{
  CharacterVector snps_from_ss;
  if(!Rf_isNull(SS) && !Rf_isMatrix(SS)){
    Function asMatrix("as.matrix");
    SS = asMatrix(SS);
  }
  if(!Rf_isNull(SS) ){
    if(!Rf_isNull(rownames(SS))){
      snps_from_ss = rownames(SS);
      vector<string> ss_names_array;
      for(int i = 0; i < snps_from_ss.length(); i++){
        ss_names_array.push_back(std::string(snps_from_ss[i]));
      }
    }
  }
  if(!Rf_isNull(SS) && X_snps.size() > 0){
    NumericMatrix SS_Matrix(SS);
    mat lp_summaries(SS_Matrix.begin(), SS_Matrix.nrow(), SS_Matrix.ncol(), false);
    Col<uword> sindex;
    if(snps_from_ss.size() > 0 ){
      snps_overlap(X_snps, snps_from_ss, xindex, sindex);
      lp_summaries = lp_summaries.rows(sindex);
    }
    if(lp_summaries.n_cols > 0 && lp_summaries.n_rows > 0){
      uvec q = find(lp_summaries < lbPval);
      lp_summaries.elem(q).fill(lbPval);
    }
    if(lp_summaries.n_rows > lp_summaries.n_cols){
      lp_summaries = lp_summaries.t();
    }
    return lp_summaries;
  }else if(!Rf_isNull(SS) && (X_snps.size() == 0 || snps_from_ss.length() == 0)){
    NumericMatrix SS_Matrix(SS);
     if(verbose){
       cout << "The SNPs' names for genotype data or Summary Statistics are "
            <<"not properly provided! " << endl;
     }
    if(P == SS_Matrix.nrow() || P == SS_Matrix.ncol()){
       if(verbose){
         cout << "The numbers of SNPs are equal for genotype data and Summary Statistics! We are going to using the p-values!" << endl;
       }
    }else{
       mat lp_summaries;
       lp_summaries.reset();
       return lp_summaries;
    }

    mat lp_summaries(SS_Matrix.begin(), SS_Matrix.nrow(), SS_Matrix.ncol(), false);
    if(lp_summaries.n_cols > 0 && lp_summaries.n_rows > 0){
      uvec q = find(lp_summaries < lbPval);
      lp_summaries.elem(q).fill(lbPval);
    }
    if(lp_summaries.n_rows > lp_summaries.n_cols){
      //to ensure each row is the summary for a GWAS
      lp_summaries = lp_summaries.t();
    }
    return lp_summaries;
  }else{
    mat lp_summaries;
    lp_summaries.reset();
    return lp_summaries;
  }

}

void convert2mat(Mat<double>*& obj, SEXP input){
  if(!Rf_isNull(input)){

    if(!Rf_isMatrix(input)){
      Function asMatrix("as.matrix");
      input = asMatrix(input);
    }

    Rcpp::NumericMatrix T(input);
    obj = new Mat<double>(T.begin(), T.nrow(), T.ncol(), false);
  }
}

void get_N_P_type(SEXP Xs, int& N, int& P, int& type){
  if(!Rf_isNull(Xs) && !Rf_isMatrix(Xs)){
    Function asMatrix("as.matrix");
    Xs = asMatrix(Xs);
  }
  if(!Rf_isInteger(Xs)){
    NumericMatrix nm(Xs);
    N = nm.nrow();
    P = nm.ncol();
    type = type_double;
  }else{
    IntegerMatrix nm(Xs);
    type = type_char;
    N = nm.nrow();
    P = nm.ncol();
  }
}

RcppExport SEXP wrap_fit(FitObj* fit){
  Rcpp::List ret;
  if(fit -> method == IGESS_METHOD){
    ret["method"] = "IGESS";
  }else if(fit -> method == LEP_METHOD){
     ret["method"] = "LEP";
  }
  ret["sigma2beta"] = fit -> sigma2beta;
  ret["sigma2e"] = fit -> sigma2e;
  ret["gammas"] = fit -> gammas;
  ret["mu"] = fit -> mu;
  ret["S"] = fit -> S;
  ret["pi"] = fit -> Pi;
  ret["M"] = fit ->P;
  ret["cov"] = fit -> cov;
  ret["L"] = fit -> L;
  ret["iter"] = fit -> iter;
  ret["u"] = fit -> u;
  ret["v"] = fit -> v;
  ret["fdr"] = 1 - fit -> gammas;
  ret["xindex"] = fit -> xindex;
  ret["time_iter"] = fit -> time_iter;
  if(fit -> pParam != NULL){
    ret["param_beta"] = (*fit -> pParam);
  }
  return ret;
}


int snps_overlap(vector<SNP>& x_snps, CharacterVector& ss_snps, Col<uword>& xindex, Col<uword>& sindex){
   vector<SNP> SS_SNPS;
   for(int i = 0; i < ss_snps.length(); i++){
     SNP snp(std::string(ss_snps[i]), i, 1);
     SS_SNPS.push_back(snp);
   }
   return snps_overlap(x_snps, SS_SNPS, xindex, sindex);
}


void preprocess_summary(SEXP Xs, SEXP  SS, double lbPval, int& type, int& N, int& P, mat& summary, void* & ptr, Col<uword>& xindex, bool verbose){
  if(!Rf_isNull(Xs) && !Rf_isMatrix(Xs)){
    Function asMatrix("as.matrix");
    Xs = asMatrix(Xs);
  }
  CharacterVector geno_snps = colnames(Xs);
  vector<SNP> xsnps;
  for(int i = 0; i < geno_snps.length(); i++){
    SNP snp(std::string(geno_snps[i]), i, 1);
    xsnps.push_back(snp);
  }
  get_N_P_type(Xs, N, P, type);
  summary = combine_summary_X(xsnps, SS, lbPval, verbose, P, xindex);
}




#endif /* lep_hpp */
