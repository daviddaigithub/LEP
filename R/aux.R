SingleVariableScan <- function(X, y){
  n <- nrow(X)
  p <- ncol(X)
  if(n != length(y)){
    cat("Not the same size of X and y")
  }
  seList = as.matrix(numeric(p));
  betaList = as.matrix(numeric(p));
  pval = as.matrix(numeric(p));
  Z = as.matrix(numeric(p));
  onesN1 = matrix(1,nrow=n,ncol=1);
  for(i in 1:p){
    x = cbind(X[,i],1);
    beta <- solve(crossprod(x),c(t(y) %*% x))
    r = y - x %*% beta
    sigma2 <- crossprod(r)/(n-1-1)
    invxtx <- solve(crossprod(x))
    Beta_cov <- c(sigma2) * invxtx
    se <- sqrt(diag(Beta_cov))
    seList[i] <- se[1]
    betaList[i] <- beta[1]
    tval <- beta / sqrt(diag(Beta_cov))
    pval_i <- 2*(1-pt(abs(tval),n - 1 -1))
    pval[i] <- pval_i[1]
    Z[i] <- tval[1]
  }
  pval[pval >= 1] = 0.99
  pval[pval <= 1e-12] = 1e-12
  return(list(betaList = betaList,seList = seList,pval = pval, Z = Z))
}


run_user_method <- function(real_model, input_method, input_summary, parallel = 10){
  draw <- real_model@draw
  pvalue <- draw@pvalue[, input_summary]
  draw@pvalue <- as.matrix(pvalue)
  innnercv <- function(i){
    sd <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    set.seed(sd)
    L <- cv_function(real_model, draw, input_method)
    result <- list(disease = real_model@name, method = input_method@name, summary = paste0(input_summary,collapse = "_"),
                   measure = "auc", value = L$auc, seed = sd, file = "filelink")
    return(as.data.frame(result))
  }

  no_cores <- parallel
  cl <- makeCluster(no_cores,
                    outfile=paste0(real_model@name,"-",
                        paste0(input_summary,collapse = "_"),"-", input_method@name,".txt"), type="FORK")
  i = 1:parallel
  results <- parLapply(cl, i, innnercv) # lapply的并行版本
  groupResult <- do.call('rbind',results) # 整合结果


  real_model@result <- c(real_model@result, new("Result", time = Sys.time(), result = groupResult))
  stopCluster(cl)
  return(groupResult)
}


copy_modelfiles_from <- function(sim, dir){
  dir <- "/home/daviddai/HKBU/LEP_Experiments/AlongU_Extend_rho"
    for(i in 1:length(sim@model_refs)){
       model <- sim@model_refs[[i]]
      # cp --parents -rf geno/files/slm/geno_1/h_0.5/h0_0.5/K_1/n_2000/n0_8000/p_10000/rho_0.8/s_0.05/u_0.05 AlongU_Extend_rho/
         system(paste0("cp "," --parents  "," -rf ", "files/",model@name, " ", dir))
    }
}


