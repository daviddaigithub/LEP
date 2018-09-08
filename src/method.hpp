//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef lep_hpp
#define lep_hpp
#include <stdio.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include "aux.hpp"
// #include <omp.h>

template <typename T>
double dotX (T* x, double* y, int n) {
    double sum = 0;
    //  #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++)
        sum += x[i] * y[i];
    return sum;
}



template<class T>
void addX (double* y, double a, T* x, int n) {
    //   #pragma omp parallel for num_threads(4)
    for (int i = 0; i < n; i++)
        y[i] += a * x[i];
}

struct PairCORAUC{
  double cor;
  double auc;
};

template<typename T>
T sum(T a1, T b1){
    return a1 + b1;
}

//template<typename T>
//LEPfit* lep(T* lpfX, vec y, int P, mat* lpsummaryinfo = NULL, Col<uword>* xindex = NULL, Options* opt = NULL, int type = type_int, bool verbose = false);

//template<typename T>
//void lep_update(T* x_j, double* gamma, double* mu, double d, double s, double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt, double* u, double* v, double* extra_vector, double* lpsummay = NULL, vec* lpparam = NULL, double xmean_j = 0, double xsd_j = 1, int type = type_int, bool befloat = true);
//

template<typename T>
void vi_update(T* x_j, double* gamma, double* mu, double d, double s, double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt, double* u, double* v, double* extra_vector, double* lpsummay = NULL, vec* lpparam = NULL, double xmean_j = 0, double xsd_j = 1, int method = IGESS_METHOD){//
    double r = (*gamma) * (*mu);
    double numerator = xy + d * r;
    double xjyhat = dotX(x_j, ytilde_pt, N);
    mat y_tilde(ytilde_pt, N, 1, false);

    bool char_or_int = strcmp(typeid(T).name(),"i") == 0 || strcmp(typeid(T).name(),"c") == 0;

    if(char_or_int)
    {
        double sum_y_tilde = as_scalar(sum(y_tilde));
        xjyhat -= xmean_j * sum_y_tilde;
        xjyhat /= xsd_j;
    }
    numerator -= xjyhat;
    (*mu) = s / sigma2e * numerator;

    double SSR_logratio = (*mu) * (*mu) / (2 * s) + 0.5 * log( s / sigma2beta);

    double ai = 0; //additional information provided by the summary statistisc
    if (lpsummay != NULL && lpparam != NULL) {
        uword K = lpparam -> size();
        double* lppara = lpparam -> memptr();
        if(method == LEP_METHOD){
          for (uword i = 0; i < K; i++) {
            double beta_f = lppara[i] * pow(lpsummay[i], lppara[i] - 1);
            ai += log((beta_f * u[i] + 1 - u[i])/(beta_f * (1 - v[i]) + v[i]));//log(lppara[i]) + (lppara[i] - 1) * log(lpsummay[i]);
          }
        }else if(method == IGESS_METHOD){
          for (uword i = 0; i < K; i++) {
            ai += log(lppara[i]) + (lppara[i] - 1) * log(lpsummay[i]);
          }
        }

        SSR_logratio += ai;
    }


    (*gamma) = 1 /(1 + exp(-(logPi + SSR_logratio)));

    double r_diff = (*gamma) * (*mu) - r;

    for(int i = 0; i < N; i++){
        extra_vector[i] = r_diff * x_j[i];
    }

    if(char_or_int)
    {
        xmean_j /= xsd_j;
        r_diff *= xmean_j;
        for (int i = 0; i < N; i++) {
            extra_vector[i] /= xsd_j;
            extra_vector[i] -= r_diff;
        }
    }
    for(int i = 0; i < N; i++){
        ytilde_pt[i] += extra_vector[i];
    }
}


template<typename T>
FitObj* fit(T* lpfX, vec y, unsigned long long P, mat* lpsummaryinfo, Col<uword>* xindex, Options* opt, int type, bool verbose, int method = IGESS_METHOD)
{
    unsigned long long N = y.n_rows;
    double mean_y;// = as_scalar(mean(y)); //mean of y, return by function 'centering'
    unsigned long long Ps = (xindex == NULL || xindex -> size() == 0) ? P : xindex -> n_elem;//P for summary
    double* mean_x = new double[Ps];
    double* sd_x = new double[Ps];
    clock_t t1 = clock();
    if(verbose){
        cout << "Start scaling the genotype matrix ...";
    }
    centering(lpfX, y, mean_y, N, P, mean_x, sd_x, xindex);
    if(verbose)
        cout <<", total time is " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;

    long long K = lpsummaryinfo != NULL ? (lpsummaryinfo -> n_rows) : 0;

    vec u;
    vec v;
    double* q11 = NULL;
    double* q00 = NULL;

    if(method == LEP_METHOD){
      u.resize(K);
      v.resize(K);
      u.fill(0.1);
      v.fill(0.9);
      q11 = new double[Ps * K];
      q00 = new double[Ps * K];
      for(int i = 0; i < Ps*K; i++){
        q11[i] = 0.1;
        q00[i] = 0.9;
      }
    }

    opt = opt != NULL ? opt : new Options();
    int max_iter = opt -> max_iter;
    int display_gap = opt -> display_gap;
    mat xty(1, Ps, fill::zeros);
    mat SZX(mean_x, 1, Ps, false);
    mat SDX(sd_x, 1, Ps, false);

    t1 = clock();

    if(verbose)
        cout << "Start calculating X'Y and diag(X'X)......";
    vec diagXTX = cal_diagXTX(y, lpfX, mean_x, sd_x, N, P, Ps, xty.memptr(), xindex);
    if(verbose)
        cout <<", time is " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;

    double pi_p = 0.01; //pi for prior proportion
    double mu0 = 0;
    double alpha0 = 0.5; //initial parameters of beta distribtuion

    Vardist vardist(Ps, mu0, pi_p);

    double sigma2e = var(y) / 2;
    double sigma2beta = sigma2e;

    vec beta = vardist.gamma % vardist.mu;
    vec ytilde = vec(N);//(*lpfX) * beta;
    ytilde.zeros();

    Col<double>* lpparams = NULL;
    if ( K > 0 ){
        lpparams = new Col<double>(K);
        lpparams -> fill(alpha0); //parameters for beta distribtuions
    }

    double L0 = -INFINITY;
    double L = 0;
    double* lpgamma = vardist.gamma.memptr();
    double* lpmu = vardist.mu.memptr();
    double* lpd = diagXTX.memptr();
    double* lpytilde = ytilde.memptr();
    double* lpxy = xty.memptr();
    double* lpsummary = lpsummaryinfo != NULL ? lpsummaryinfo -> memptr() : NULL;
    uword iter = 0;
    mat extra_vector(N, 1, fill::zeros);
    clock_t t_iteration = clock();
    for (iter = 0; iter < max_iter; iter ++) {
        clock_t t1 = clock();
        if (verbose) {
            if(iter == 0)  cout <<"Begin Iterations" << endl;
        }
        double logPi = log(pi_p / (1 - pi_p));
        double sigma2e_Sigma2beta = sigma2e / sigma2beta;
        vec xxsigma = diagXTX + sigma2e_Sigma2beta;
        vardist.sigma2beta = sigma2e / xxsigma;
        double* S = vardist.sigma2beta.memptr();
        double gamma_sum = 0;
        if(xindex == NULL || xindex -> size() == 0){
            for (uword j = 0; j < P; j++) {
                vi_update(lpfX + N * j, lpgamma + j, lpmu + j, lpd[j], S[j],
                          logPi, sigma2beta, sigma2e, N,  lpxy[j], lpytilde,
                            u.memptr(), v.memptr(), extra_vector.memptr(),
                               lpsummary + K * j, lpparams, mean_x[j], sd_x[j], method);
                gamma_sum += *(lpgamma + j);
            }
        }else{
            uword j = 0;
            for (uword i = 0; i < xindex -> size(); i++){
                j = xindex -> at(i);
                vi_update(lpfX + N * j, lpgamma + i, lpmu + i, lpd[i], S[i],
                          logPi, sigma2beta, sigma2e, (int)N,  lpxy[i], lpytilde,
                          u.memptr(), v.memptr(), extra_vector.memptr(),
                                   lpsummary + K * i, lpparams, mean_x[i], sd_x[i], method);
                gamma_sum += *(lpgamma + i);
            }

        }
        //update alpha_k for beta-distribution
        update_betaparam(Ps, K, gamma_sum, lpsummary, lpgamma, lpparams, q11, q00, xindex, method);


        // update sigma2beta, sigma2e, u, v
        update_param(N, Ps, vardist, sum(square(y-ytilde)), diagXTX, sigma2e, sigma2beta, pi_p, q11, q00, u.memptr(), v.memptr(), (int)K, lpgamma, gamma_sum, xindex, method);


        L = lb_linear(ytilde, diagXTX,y, sigma2e, vardist) + lb_gamma(vardist.gamma, logPi) +lb_klbeta(vardist, sigma2beta)
        + lb_pvalue(Ps, K, lpsummary, lpgamma, lpparams, q11, q00, u.memptr(), v.memptr(), method);

        if(method == LEP_METHOD){
          update_q11_q00(Ps, K, lpsummary, lpparams, u.memptr(), v.memptr(), q11, q00);
        }

        if(verbose && (iter % display_gap == 0)){
            cout <<iter <<"th iteration L=" << L << ";sigma2e = " << sigma2e << "; sigma2beta = " << sigma2beta << " ;pi = " <<
            pi_p <<" time = " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
            cout <<"uv for " << K <<" GWAS" << endl;
            cout <<"u =" << u.t() << endl;
            cout <<"v =" << v.t() << endl;
        }

        if(L < L0){
            cout << "Lowerbound decreasing, Error at iteration "<<iter <<"th iteration, diff = " << L-L0 << endl;
            break;
        }else if(L - L0 < 1e-5)//if(fabs((L - L0)/L0) < 1e-8) //
        {
            if(verbose){
                cout<<"Converge at " << iter << "th iteration, L = "<< L << endl;
            }
            break;
        }
        L0 = L;

    }
    double time_iter = (clock() - t_iteration) * 1.0 / CLOCKS_PER_SEC;
    if(verbose){
      cout << "The total time for " << iter << "  iterations is " << time_iter << endl;
    }


    mat sd_vec(sd_x, Ps, 1, false);
    vardist.mu /= sd_vec;
    mat betahat = conv_to<mat>::from(vardist.gamma % vardist.mu);
    double cov = mean_y - as_scalar(SZX * betahat);
    FitObj* fit = new FitObj(N, Ps,  K, iter, L,  sigma2e, sigma2beta, pi_p, vardist.gamma, vardist.mu
                             , vardist.sigma2beta, lpparams,  ytilde, cov, u, v, xindex, method);
    fit -> time_iter = time_iter;
    delete[] mean_x;
    delete[] sd_x;
    return fit;
}
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold);


template<class T>
double CV_Proc(T* lpfX, vec y, int P, mat* lpsummaryinfo = NULL, Col<uword>* xindex = NULL, Options* opt = NULL,int type = type_int, std::string measure = "auc", bool verbose = true, int method = IGESS_METHOD)
{
  uword N = y.n_rows;
  opt = opt != NULL ? opt : new Options();
  uword nfold = opt -> n_fold;
  Col<uword> indices = cross_valind(N, nfold);
  vec ylabel = y;
  vec predY(N);
  if (verbose) {
    cout << "Start " << nfold <<" cross validations!" << endl;
  }
  bool befloat = true;
  for (uword i = 1; i <= nfold; i++) {
    if (verbose) cout <<".";
    Col<uword> train_idx = find(indices != i);
    Col<uword> test_idx = find(indices == i);
    Mat<T> mat_f(lpfX, N, P);
    Mat<T> trainM = mat_f.rows(train_idx);
    vec ytrain = y(train_idx);
    Mat<T> testM = mat_f.rows(test_idx);
    vec ytest = y(test_idx);
    T* trainX = trainM.memptr();
    FitObj* f = fit(trainX, ytrain, P,  lpsummaryinfo, xindex, opt, type, verbose, method);
    vec predy = f -> predict(&testM);
    predY.elem(test_idx) = predy;
    delete f;
  }

  double predict = 0;
  if(measure.compare("auc") == 0 ){
    predict = calauc(conv_to<vec>::from(ylabel), conv_to<vec>::from(predY));
  }else if(measure.compare("mse") == 0 ){
    vec diff = y - predY;
    predict = mean(diff % diff);
  }else if(measure.compare("cor") == 0 ){
    predict = as_scalar(cor(y, predY));
  }
  return predict;
}






//template<class T>
//void addX (double* y, double a, T* x, int n);
//
//template <typename T>
//double dotX (T* x, double* y, int n);
#endif /* lep_hpp */
