////  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef aux_hpp
#define aux_hpp

#include <stdio.h>
#include <armadillo>

#include "readPlink.hpp"
using namespace std;
using namespace arma;
//#include <omp.h>

#define type_int 2
#define type_char 0
#define type_double 1

#define LEP_METHOD 0
#define IGESS_METHOD 1

const int kPValueVer = 0;
const int kVbVer = 1;

class PerformanceObject{//object to record the performance
public:
    double FDR;
    double power;
    double auc;

};

class Vardist{
public:
    vec gamma;
    vec mu;
    vec sigma2beta;
    Vardist(uword P, double mu0, double gamma0){
        gamma = vec(P);
        gamma.fill(gamma0);
        mu = vec(P);
        mu.fill(0);
    }
};

class DataModel{ //model to store the original data
public:
    DataModel(mat* X, vec* y, Col<int>* labels, vec* beta, double h){
        this -> X = new Mat<double>(X -> memptr(), X->n_rows, X -> n_cols, false);
        this -> y = new Col<double>(y -> memptr(), y->size(), false);
        this -> labels = new Col<int>(labels -> memptr(), labels->size(), false);
        this -> beta = new Col<double> (beta -> memptr(), beta -> size(), false);
        this -> h = h;

    };
    ~DataModel(){
        delete this -> X;
        delete this -> y;
        delete this -> labels;
        delete this -> beta;
    }
    mat* X;
    vec* y;
    Col<int>* labels;
    vec* beta;
    double h; // heritability

};


class FitObj{
    //class for the model generated

public:
  FitObj( uword N, uword P,  uword K, uword iter, double L,  double sigam2e, double sigma2beta, double Pi, vec gammas, vec mu, vec S, vec* pParam, vec Xr, double cov, vec u_, vec v_, Col<uword>* xindex, int method = IGESS_METHOD){
        this -> N = N;
        this -> P = P;
        this -> K = K;
        this -> L = L;
        this -> iter = iter;

        this -> sigma2e = sigam2e;
        this -> sigma2beta = sigma2beta;
        this -> Pi = Pi;

        this -> gammas = gammas;
        this -> mu = mu;
        this -> S = S;
        this -> pParam = pParam;
        this -> Xr = Xr;
        this -> cov = cov;
        this -> u = u_;
        this -> v = v_;
        if(xindex == NULL || xindex -> size() == 0){
          Col<uword> xindex0;
          this -> xindex = xindex0;
        }else{
          this -> xindex = *xindex;
        }
        this -> method = method;
        // this -> gamma_K = mat(gamma_K, P, K);

    }

  FitObj(vec gammas, vec mu, double cov){
       this -> gammas = gammas;
       this -> mu = mu;
       this -> cov = cov;
    }
    ~FitObj( ){

    }

    // vec predict( mat* X );
    template<class T>
    vec predict(Mat<T>* X)
    {
      vec beta;
      vec yhat;
      if(this -> xindex.size() == 0){
        beta = conv_to<vec>::from(this->gammas % this->mu);
      }else if(this -> xindex.size() < X ->n_cols){
        beta.resize(X ->n_cols);
        beta.zeros();
        beta(this -> xindex) = conv_to<vec>::from(this->gammas % this->mu);
      }
      yhat = this->cov + (*X) * beta;
      return yhat;
    }

    void cal_powerfdr(DataModel* model, double threshold,PerformanceObject* po);
    double cal_auc(DataModel* model);

    double sigma2e;
    double sigma2beta;
    double Pi; //number of variables
    double L; //lower bound
    uword N;
    uword P;
    uword K;
    uword iter;
    vec gammas;
    vec mu;
    vec S;
    vec* pParam;
    vec Xr;
    Col<uword> xindex;
    double cov;
    vec u;
    vec v;
    int method;
    double time_iter;
    // mat gamma_K;
};

class Options{
public:
    Options(){
        this -> max_iter = 300;
        this -> display_gap = 60;
        this -> n_fold = 5;
    }

    Options(int max_iter, int display_gap){
        this -> max_iter = max_iter;
        this -> display_gap = display_gap;
    }

    Options(int max_iter, int display_gap, int n_fold){
        this -> max_iter = max_iter;
        this -> display_gap = display_gap;
        this -> n_fold = n_fold;
    }
    int max_iter;
    int display_gap;
    int n_fold;

};


vec fdr2FDR(vec fdr);

double lb_linear(vec ytilde, vec diagXTX, vec y, double sigma2e, Vardist vardist);
double lb_gamma(vec gamma, double log_pi);
double lb_klbeta(Vardist vardist, double sigma2beta);


//double dotX (double* x, double* y, int n);






/*update the parameters of the gamma distributions*/
// void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams);
// void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00);

// double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams);
double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams,double* q11, double* q00, double* u, double* v, int method = IGESS_METHOD);

// void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p);
void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p, double* q11, double* q00, double* u, double* v, int K, double* lpgamma, double gamma_sum, Col<uword>* xindex = NULL, int method = IGESS_METHOD);

double clog2(double v);
//remove the effects of covariates Z
template <typename T>
void calMeanSd(T * data, int N, int M, double* mean, double* sd, Col<uword>* xindex){
    int value;
    if(xindex == NULL || xindex -> size() == 0){
        for(int j = 0; j < M; j++){
            mean[j] = 0;
            sd[j] = 0;
            T* col_j = data + (long long) j * (long long)  N;
            for (int i = 0; i < N; i++) {
                value = col_j[i];//data[j * N + i];
                mean[j] += value;
                sd[j] += value * value;
            }
            sd[j] = sqrt(sd[j] / N - (mean[j] / N) * (mean[j] / N));
            mean[j] = mean[j] / N;
        }
    }else{
        uword j = 0;
        for(int k = 0; k < xindex -> size(); k++){
            j = xindex -> at(k);
            mean[k] = 0;
            sd[k] = 0;
            T* col_j = data + (long long) j * (long long)  N;
            for (int i = 0; i < N; i++) {
                value = col_j[i];//data[j * N + i];
                mean[k] += value;
                sd[k] += value * value;
            }
            sd[k] = sqrt(sd[k] / N - (mean[k] / N) * (mean[k] / N));
            mean[k] = mean[k] / N;
        }

    }

}

double calauc(vec label, vec pred);

template<class T>
void centering(T* X, vec& y, double& mean_y, long long N, long long P, double* x_mean, double* x_sd, Col<uword>* xindex){
    calMeanSd(X, N, P, x_mean, x_sd, xindex);
    mean_y = as_scalar(mean(y)); //mean of y
    y -= mean_y;
    //T is double
    if(strcmp(typeid(T).name(),"i") != 0 && strcmp(typeid(T).name(),"c") != 0){
      if(xindex == NULL || xindex -> size() == 0){
        for (uword j = 0; j < P; j++) {
            T* col_j = X + j * N;
            for(uword i = 0; i < N; i++){
                col_j[i] -= x_mean[j];
                col_j[i] /= x_sd[j];
            }
        }
      }else{
        long long j = 0;
        for(int k = 0; k < xindex -> size(); k++){
          j = (long long)xindex -> at(k);
          T* col_j = X + j * N;
          for (int i = 0; i < N; i++){
            col_j[i] -= x_mean[k];
            col_j[i] /= x_sd[k];
          }
        }
      }

    }
}

template<class T>
void ytX(double* y, T* data, long long N, long long M, double* xty, Col<uword>* xindex){
    if(xindex == NULL || xindex -> size() == 0){
        for(long long j = 0; j < M; j++){
            xty[j] = 0;
            T* col_j = data + j * N;
            for (int i = 0; i < N; i++) {
                xty[j] += y[i] * col_j[i];

            }
        }
    }else{
        long long j = 0;
        for(long long k = 0; k < xindex -> size(); k++){
            xty[k] = 0;
            j = xindex -> at(k);
            T* col_j = data + j * N;
            for (int i = 0; i < N; i++) {
                xty[k] += y[i] * col_j[i];
            }
        }
    }

}

template<class T>
vec cal_diagXTX(vec& y, T* x, double* SZX, double* SDX, unsigned long long N,unsigned long long M, unsigned long long Ms, double* xty,  Col<uword>* xindex){
    double sum_y = sum(y);
    vec diagXTX;
    ytX(y.memptr(), x, N, M, xty, xindex);
    if(strcmp(typeid(T).name(),"i") == 0 || strcmp(typeid(T).name(),"c") == 0)
    {
        for(int i = 0; i < Ms; i++){
            xty[i] -= sum_y * SZX[i];
            xty[i] /= SDX[i];
        }
    }
    mat diag(Ms,1);
    diag.fill(N);
    diagXTX = diag;
    return diagXTX;
}


/*update the parameters of the gamma distributions*/
void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00, Col<uword>* xindex, int method = IGESS_METHOD);
void update_q11_q00(uword P, uword K, double * lpsummary,vec* lpparams, double* u, double* v, double* q11, double* q00);




#endif /* lep_aux_hpp */
