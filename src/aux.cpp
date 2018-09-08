//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#include "aux.hpp"

double calauc(arma::vec label, arma::vec pred){
    double auc = 0;
    double m = mean(label);
    vec label2 = label;
    label2(find(label >= m)).fill(1);
    label2(find(label <= m)).fill(0);
    label = label2;
    uword N = pred.size();
    uword N_pos = sum(label);
    uvec  idx = sort_index(pred,"descend");
    vec above = linspace<vec>(1, N, N) - cumsum(label(idx));
    auc = (1 - sum(above % label(idx)) / (N_pos * (N-N_pos)));
    auc = auc > 0.5?auc:(1-auc);
    return auc;
}

/*update the parameters of the gamma distributions*/
void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00, Col<uword>* xindex, int method)
{
    if(K == 0) return;
    vec alphalogpvec(K);
    alphalogpvec.fill(0);
    double* lpalphalogp = alphalogpvec.memptr();

    for(uword k=0; k < K; k++)
    {
        gamma_sum = 0;

        if(xindex == NULL || xindex -> size() == 0){
           if(method == LEP_METHOD){
             for (uword j = 0; j < P; j++)
             {
               double pvalue = lpsummary[K * j + k];
               double coef = lpgamma[j] * q11[P*k + j] + (1 - lpgamma[j]) * (1 - q00[P*k + j]);
               gamma_sum += coef;
               lpalphalogp[k] += coef * ( -log(pvalue) );
             }
           }else if(method == IGESS_METHOD){
             for (uword j = 0; j < P; j++)
             {
               double pvalue = *(lpsummary + K * j + k);
               gamma_sum += lpgamma[j];
               lpalphalogp[k] += (lpgamma[j])*(-log(pvalue));
             }
           }

        }else{
          if(method == LEP_METHOD){
            uword j = 0;
            for(int i = 0; i < xindex -> size(); i++){
              j = xindex -> at(i);
              double pvalue = lpsummary[K * i + k];
              double coef = lpgamma[i] * q11[P*k + i] + (1 - lpgamma[i]) * (1 - q00[P*k + i]);
              gamma_sum += coef;
              lpalphalogp[k] += coef * ( -log(pvalue) );
            }
          }else if(method == IGESS_METHOD){
            for(int i = 0; i < xindex -> size(); i++){
              double pvalue = lpsummary[K * i + k];
              gamma_sum += lpgamma[i];
              lpalphalogp[k] += lpgamma[i] * (-log(pvalue));
            }
          }


        }
        (*lpparams)[k] = gamma_sum / lpalphalogp[k];
        if( (*lpparams)[k] > 1){
          (*lpparams)[k] = 1;
        }
    }

}




void update_q11_q00(uword P, uword K, double * lpsummary,vec* lpparams, double* u, double* v, double* q11, double* q00)
{
  if(K==0) return;
  double* lpparam = lpparams -> memptr();
  for(int j = 0; j < P; j++)
    for(int k = 0; k < K; k++){
      double param = lpparam[k];
      double pvalue = lpsummary[K * j + k];
      double beta_f = param * pow(pvalue, param - 1);
      q11[P*k + j] = beta_f * u[k] / (beta_f * u[k] + 1 - u[k]);
      q00[P*k + j] = v[k] / (beta_f * (1 -v[k]) +  v[k]);
    }
}

void FitObj::cal_powerfdr(DataModel* model, double threshold, PerformanceObject* po)
{
    vec gFDR = fdr2FDR(this->gammas);
    uvec ufound = find(gFDR < threshold);
    uword nerr = sum((*model -> labels)(ufound) == 0);
    uword nfound = ufound.size();

    double FDR = (nfound != 0) ? nerr * 1.0 / nfound : 0;

    uword nCausal = sum((*model ->labels) != 0);
    double power =  (ufound.size() - nerr) * 1.0  / nCausal ;

    po -> FDR = FDR;
    po -> power = power;

}

double FitObj::cal_auc(DataModel* model){
    vec label = conv_to<vec>::from((*model -> labels));
    double auc = calauc(label, this -> gammas);
    return auc;
}


double clog2(double v){
  return log(v + (v < 1e-12));
}

double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00, double* u, double* v, int method){
  double lb = 0;
  if(K == 0) return lb;
  double* lpparam = lpparams -> memptr();
  double pvalue = 0;
  double param_k = 0;
  if(method == IGESS_METHOD){
    for (int j = 0; j < P; j++) {
      for (int k = 0; k < K; k++) {
        param_k = lpparam[k];
        pvalue =  lpsummary[j * K + k];
        lb += lpgamma[j] * (log(param_k) + (param_k - 1)*log(pvalue));
      }
    }
  }else if(method == LEP_METHOD){
    for (int j = 0; j < P; j++) {
      double sum_q1 = 0;
      double sum_q0 = 0;
      for (int k = 0; k < K; k++) {
        param_k = lpparam[k];
        pvalue =  lpsummary[j * K + k];
        double q11_jk = q11[P*k + j];
        double q10_jk = 1 - q11_jk;
        double q00_jk = q00[P*k + j];
        double q01_jk = 1 - q00_jk;

        sum_q1 += q11_jk * (clog2(param_k) +  (param_k -1)*clog2(pvalue) + clog2(u[k]) - clog2(q11_jk)) + q10_jk * (clog2(1 - u[k]) - clog2(q10_jk));
        sum_q0 += q01_jk * (clog2(param_k) +  (param_k -1)*clog2(pvalue) + clog2(1 -v[k]) - clog2 (q01_jk)) ;
        sum_q0 += q00_jk * (clog2(v[k])- clog2(q00_jk));

      }
      lb +=  lpgamma[j] * sum_q1 + (1 - lpgamma[j]) * sum_q0;
    }
  }

  return lb;
}

void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p, double* q11, double* q00, double* u, double* v, int K, double* lpgamma, double gamma_sum, Col<uword>* xindex, int method){
    vec term1 = vardist.gamma % (vardist.sigma2beta + square(vardist.mu));
    double term2 = sum((term1 - square(vardist.gamma % vardist.mu)) % diagXTX);

    sigma2e = (sumyytilde + term2) / N;
    double sum_vardist_gamma = sum(vardist.gamma);
    pi_p = sum_vardist_gamma / P;
    sigma2beta = sum(term1) / sum_vardist_gamma;

    if(method == IGESS_METHOD) return;

    for(int k = 0; k < K; k++){

        double sum_q11 = 0;
        double sum_q00 = 0;

        if(xindex == NULL || xindex -> size() == 0){
            for(int j = 0; j < P; j++){
                sum_q11 += lpgamma[j] * q11[P * k + j];
                sum_q00 += (1 - lpgamma[j]) * q00[P * k + j];
            }
        }else{
            uword j = 0;
            for(int i = 0; i < xindex -> size(); i++){
                sum_q11 += lpgamma[i] * q11[P*k + i];
                sum_q00 += (1 - lpgamma[i]) * q00[P*k + i];
            }
        }

        u[k] = sum_q11 / gamma_sum;
        v[k] = sum_q00 / (P - gamma_sum);
    }
}

double dotXX (double* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}

double dotXX (float* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}


//lower bound for the linear part
double lb_linear(vec ytilde, vec diagXTX, vec y, double sigma2e, Vardist vardist){
    uword n = y.n_elem;
    double lb1 = - (0.5*n)*log(2 * M_PI * sigma2e);
    double lb2 = - 0.5 * sum(square(y-ytilde))/sigma2e;
    double lb3 =
    -0.5*sum( (vardist.gamma % (vardist.sigma2beta + square(vardist.mu)) -
               square(vardist.gamma % vardist.mu)) % diagXTX )/sigma2e;
    return lb1 + lb2 + lb3;

};

vec logpexp(vec x) {
    vec y = x;
    uvec idx = (x < 16);
    y(idx) = log(1 + exp(x(idx)));
    return y;
}

double logpexp(double x) {
    double y = log(1 + exp(x));
    return y;
}

//lower bound for gamma part
double lb_gamma(vec gamma, double log_pi){
    return sum((gamma - 1) * log_pi + (-logpexp(-log_pi)));
};

//lower bound for the klbeta part
double lb_klbeta(Vardist vardist, double sigma2beta){
    double lb = -sum(vardist.gamma % log(vardist.gamma+(vardist.gamma==0)) + (1-vardist.gamma) % log(1-vardist.gamma+(vardist.gamma == 1))) + 0.5*sum(vardist.gamma % (1+log(vardist.sigma2beta / sigma2beta)- (square(vardist.mu) + vardist.sigma2beta)/ sigma2beta ));
    return lb;
}

//convert local fdr to Global FDR
vec fdr2FDR(vec fdr){
    uword M = fdr.size();
    uvec indices = sort_index(fdr);
    vec sort_fdr = fdr(indices);
    vec FDR = cumsum(sort_fdr) / linspace<vec>(1, M, M);
    FDR.elem(find(FDR  > 1)).ones();
    FDR(indices) = FDR;
    return FDR;
}
















