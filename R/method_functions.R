Lasso <- new_method("Lasso", "Lasso Path Algorithm", method = function(model, draw){
  if("cc" %in% names(model@params)){
    x <- model@params$x[draw@train_index,]
    tx <- model@params$x[draw@test_index,]
  }else{
    x <- model@params$x
    tx <- model@params$tx
  }
  cvfit <- cv.glmnet(x, draw@y)
  tyhat <- predict(cvfit, tx, s = "lambda.min")
  coef0 <- coef(cvfit, s = "lambda.min")
  beta <- coef0[2:nrow(coef0)]
  list(beta = beta, tyhat = tyhat)
})

LassoCC <- new_method("LassoB", "Lasso Binomial Family", method = function(model, draw){
  if("cc" %in% names(model@params)){
    x <- model@params$x[draw@train_index,]
    tx <- model@params$x[draw@test_index,]
  }else{
    x <- model@params$x
    tx <- model@params$tx
  }
  cvfit <- cv.glmnet(x, draw@y,family="binomial")
  tyhat <- predict(cvfit, tx, s = "lambda.min")
  coef0 <- coef(cvfit, s = "lambda.min")
  beta <- coef0[2:nrow(coef0)]
  list(beta = beta, tyhat = tyhat)
})



make_folds <- function(n, nfolds) {
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}

bvsr <- new_method("BVSR", "Bayesian Variable Selection Regression",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      obj <- LEP::LEP(x, draw@y, opts=list(max_iter=300,display_gap=1,n_fold=5), verbose = F)
                      list(beta = obj$mu * obj$gamma,
                           tyhat = tx %*% (obj$ mu * obj$gamma) + as.numeric(obj$cov),fdr = obj$fdr, cov = as.numeric(obj$cov),fit = obj)
                    })

igess <- new_method("IGESS", "IGESS Algorithm",
                     method = function(model, draw) {
                       if("cc" %in% names(model@params)){
                         x <- model@params$x[draw@train_index,]
                         tx <- model@params$x[draw@test_index,]
                       }else{
                         x <- model@params$x
                         tx <- model@params$tx
                       }
                       obj <- IGESS::IGESS(x, draw@y, SS = draw@pvalue,opts=list(max_iter=300,display_gap=1,n_fold=5))
                       list(beta = obj$mu * obj$gamma,
                            tyhat = tx %*% (obj$ mu * obj$gamma) + as.numeric(obj$cov),
                            fdr = obj$fdr, cov = as.numeric(obj$cov))
                     })


lep <- new_method("LEP", "LEP Algorithm",
                    method = function(model, draw) {
                     if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                     }else{
                        x <- model@params$x
                        tx <- model@params$tx
                     }
                     obj <- LEP::LEP(x, draw@y, SS = draw@pvalue, opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                     result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                     return(result)
                    })
lep100 <- new_method("LEP100", "LEP Algorithm taking the pvalues of the first column",
                   method = function(model, draw) {
                     if("cc" %in% names(model@params)){
                       x <- model@params$x[draw@train_index,]
                       tx <- model@params$x[draw@test_index,]
                     }else{
                       x <- model@params$x
                       tx <- model@params$tx
                     }
                     PV <- draw@pvalue[,1]
                     obj <- LEP::LEP(x, draw@y, SS = cbind(PV, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                     result <- list(beta =  (obj$ mu * obj$gammas),
                                    tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                     return(result)
                   })

pval_permute <- function(pval, ratio){
  if(ratio < 0 || ratio > 1){
    ratio <- 1;
  }
  PV <- pval
  P <- length(PV)
  I <- sample.int(P, size = round(ratio*P), replace = FALSE);
  PV0 <- PV[I]
  PV0 <- sample(PV0, round(ratio*P), replace = FALSE);
  PV[I] <- PV0
  return(PV)
}

lep90 <- new_method("LEP90", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,0.1);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })

lep70 <- new_method("LEP70", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,0.3);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })

lep50 <- new_method("LEP50", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,0.5);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })

lep30 <- new_method("LEP30", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,0.7);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })

lep10 <- new_method("LEP10", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,0.9);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })

lep00 <- new_method("LEP00", "LEP Algorithm taking the pvalues of the first column",
                    method = function(model, draw) {
                      if("cc" %in% names(model@params)){
                        x <- model@params$x[draw@train_index,]
                        tx <- model@params$x[draw@test_index,]
                      }else{
                        x <- model@params$x
                        tx <- model@params$tx
                      }
                      pval <- draw@pvalue[,1]
                      PV <- pval_permute(pval,1);
                      obj <- LEP::LEP(x, draw@y, SS = cbind(pval, PV), opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                      result <- list(beta =  (obj$ mu * obj$gammas),
                                     tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                      return(result)
                    })


lep0 <- new_method("LEP0", "LEP Algorithm taking the pvalues of the first column",
                  method = function(model, draw) {
                    if("cc" %in% names(model@params)){
                      x <- model@params$x[draw@train_index,]
                      tx <- model@params$x[draw@test_index,]
                    }else{
                      x <- model@params$x
                      tx <- model@params$tx
                    }
                    obj <- LEP::LEP(x, draw@y, SS = draw@pvalue[,1], opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                    result <- list(beta =  (obj$ mu * obj$gammas),
                                   tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                    return(result)
                  })

lep1 <- new_method("LEP1", "LEP Algorithm taking the pvalues of the relevant columns",
                   method = function(model, draw) {
                     if("cc" %in% names(model@params)){
                       x <- model@params$x[draw@train_index,]
                       tx <- model@params$x[draw@test_index,]
                     }else{
                       x <- model@params$x
                       tx <- model@params$tx
                     }
                     K = model@params$K;
                     obj <- LEP::LEP(x, draw@y, SS = draw@pvalue[,1:K], opts=list(max_iter=300,display_gap=1,n_fold=5), verbose =  F)
                     result <- list(beta =  (obj$ mu * obj$gammas),
                                    tyhat = tx %*% (obj$ mu * obj$gammas) + as.numeric(obj$cov), cov = as.numeric(obj$cov), fdr = obj$fdr, u = obj$u, v = obj$v,fit = obj)
                     return(result)
                   })




gpa <- new_method("GPA", "Genetic Analysis Incorporating Pleiotropy and Annotation",
                  method = function(model, draw) {
                    if("cc" %in% names(model@params)){
                      x <- model@params$x[draw@train_index,]
                      tx <- model@params$x[draw@test_index,]
                    }else{
                      x <- model@params$x
                      tx <- model@params$tx
                    }
                    n = nrow(x)
                    p = ncol(x)
                    obj <- GPA::GPA(cbind(draw@pvalue0, draw@pvalue[,1]))
                    list(tyhat = matrix(0,n,1), fdr = obj@fit$Z[,1] + obj@fit$Z[,3], beta = matrix(0, p, 1))
                  })

pvalue <- new_method("pvalue", "pvalue",
                     method = function(model, draw) {
                       obj <- GPA::GPA(cbind(draw@pvalue0))
                       # save(obj,file="GPAObj.RData")
                       list(fdr = obj@fit$Z[,1])
                     })




cv <- new_method_extension("cv", "cross validated",
                           method_extension = function(model, draw, out,
                                                       base_method) {
                             nfolds <- 5
                             Y <- matrix(NA, model$n, 1)
                             ii <- make_folds(model$n, nfolds)
                             for (i in seq_along(ii)) {
                               train <- model
                               train@params$x <- model@params$x[-ii[[i]], ]
                               train@params$n <- nrow(train@params$x)

                               train_draw <- draw
                               train_draw@y = as.matrix(draw@y[-ii[[i]],])

                               test <- model
                               test@params$x <- model@params$x[ii[[i]], ]
                               test@params$n <- nrow(test@params$x)
                               fit <- base_method@method(model = train,
                                                         draw = train_draw)
                               yhat <- test@params$x %*% fit$beta + fit$cov
                               Y[ii[[i]]] <- yhat

                             }
                             list(cvy = Y)
                           })

