library(simulator)

makcov2 <- function(p, rho){
  sigma <- matrix(0,p,p)
  sigma <- rho^abs(row(sigma)-col(sigma))
  return(sigma)
}

simuGamma <- function(P, K, rho){
  sigma <- makcov2(K, rho)
  L <- chol(sigma)
  sparsity <- 0.05
  Gammas <- matrix(rnorm(P * K), P, K)
  Gammas <- Gammas %*% L
  maf <- matrix(qnorm(sparsity), P, K)
  G_out <- matrix(0, P, K)
  G_out[Gammas < maf] = 1
  return(G_out)
}

simuGamma2 <- function(P, K, rho, sigma){
  # sigma <- makcov2(K, rho)
  L <- chol(sigma)
  sparsity <- 0.05
  Gammas <- matrix(rnorm(P * K), P, K)
  Gammas <- Gammas %*% L
  maf <- matrix(qnorm(sparsity), P, K)
  G_out <- matrix(0, P, K)
  G_out[Gammas < maf] = 1
  return(G_out)
}

getbeta <- function(x, p, h, s, se){
  nc <- round(p * s) #number of associated SNPs
  assoc = sample.int(p, nc)
  beta <- matrix(0,p,1)
  b <- sqrt(1e-2) * rnorm(nc);
  beta[assoc] <- b
  sb = h/(1-h)/var(x %*% beta);
  beta = as.numeric(sqrt(sb*se))*beta;
  return(beta)
}

getbeta_give_gamma <- function(x, p, h, gamma, se){
  nc <- sum(gamma != 0) #number of associated SNPs
  beta <- matrix(0,p,1)
  b <- sqrt(1e-2) * rnorm(nc);
  beta[gamma != 0] <- b
  sb = h/(1-h)/var(x %*% beta);
  beta = as.numeric(sqrt(sb*se))*beta;
  return(beta)
}

getbeta4p <- function(x, p, h, coef, se){
  #input h, p, s, x, se
  #output beta
  beta <- matrix(0,p,1)
  b <- sqrt(1e-2) * rnorm(sum(coef == 1));
  beta[coef == 1] <- b
  sb = h/(1-h)/var(x %*% beta);
  beta = as.numeric(sqrt(sb*se)) * beta;
  return(beta)
}

generate_pvalue <- function(n, p, h, rho, coef0, u, v, geno, maf){
     if(geno == 0){
       x <- createX(n, p, rho)
     }
     else{
       # maf = 0.05 + 0.45*matrix(runif(p),1,p);
       x <- createGeno(n, p, rho, maf)
     }
     coef <- numeric(p)
     if(u > 0 && v > 0){
       coef[coef0 == 1] = (runif(sum(coef0 == 1)) < u)
       coef[coef0 == 0] = (runif(sum(coef0 == 0)) < 1 - v)
     }else{
       nz = sum(coef0 == 1);
       coef[sample.int(p, nz)] = 1
     }
     se <- 4
     beta <- getbeta4p(x, p, h, coef, se)
     y <- x %*% beta + sqrt(se) * rnorm(n);
     obj <- SingleVariableScan(x, y)
     return(list(p=obj$pval, gamma = coef))
}

setClass("Draw", representation(y = "matrix", pvalue = "matrix", pvalue0 = "matrix", ty = "matrix", train_index = "numeric", test_index = "numeric", gamma_one = "matrix"))

#generate gaussian distributed X
createX <- function(n, p, rho){
  x <- matrix(rnorm(n*p),n,p)
  if(rho == 0 || rho < 0) return(x);
  blocksize <- 1000
  sigma <- makcov2(blocksize, rho)
  nb <- p %/% blocksize
  left <- p - nb * blocksize
  l <- chol(sigma)
  L <- matrix(0, p, p)
  for(i in 1:nb){
    istart = ((i-1)*blocksize + 1)
    iend = i*blocksize
    L[istart : iend, istart : iend] = l
  }
  if(left > 0){
    istart = nb * blocksize + 1
    L[istart : p, istart : p] = l[1:left,1:left]
  }
  x <- x %*% L;
  return(x)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

createGeno <- function(n, p, rho, maf){
  blocksize <- 1000
  sigma <- makcov2(blocksize, rho)
  nb <- p %/% blocksize
  left <- p - nb * blocksize
  l <- chol(sigma)
  L <- matrix(0, p, p)
  for(i in 1:nb){
    istart = ((i-1)*blocksize + 1)
    iend = i*blocksize
    L[istart : iend, istart : iend] = l
  }
  if(left > 0){
    istart = nb * blocksize + 1
    L[istart : p, istart : p] = l[1:left,1:left]
  }

  R1 <- matrix(rnorm(n*p),n,p);
  R1 <- R1 %*% L
  R2 <- matrix(rnorm(n*p),n,p);
  R2 <- R2 %*% L
  MAFMat = rep.row(qnorm(maf),n)
  X1 = matrix(0,n,p)
  X2 = matrix(0,n,p)
  X1[R1 < MAFMat] = 1
  X2[R2 < MAFMat] = 1
  x <- X1 + X2
  return(x)
}

get_case_control_samples <- function(n, geno, p, h, rho, s, u){
  N = (n/2)*20 + 1000
  if(geno == 0){
    x <- createX(N, p, rho)
  }else{
    maf = 0.05 + 0.45*matrix(runif(p),1,p);
    x <- createGeno(N, p, rho,maf)
  }
  se <- 4
  beta <- getbeta(x, p, h, s, se = se)
  mu <- as.numeric(x %*% beta)
  return(list(x = x, mu = mu, beta = beta, N = N))
}

simulate_y <- function(mu, se, n, N){
  y <- as.matrix(mu + sqrt(se) * rnorm(N))
  stdY = (y - mean(y))/sd(y);
  cc = 0.1;
  cutZ = VGAM::probit(1-cc);
  y = as.double(stdY > cutZ)
  nonzeros <- which(y==1)
  nz <- length(nonzeros)
  nonzeros <- nonzeros[sample.int(nz,nz,replace = F)]
  zeros <- which(y==0)
  z <- length(zeros)
  zeros <- zeros[sample.int(z,z,replace = F)]

  train_pos1 <- nonzeros[1:(n/2)]
  train_pos0 <- zeros[1:(n/2)]
  test_pos1 <- nonzeros[(n/2+1):n]
  test_pos0 <- zeros[(n/2+1):n]

  train_index <- c(train_pos1,train_pos0)
  test_index <- c(test_pos1, test_pos0)
  Y <- y[train_index]
  tY <- y[test_index]
  Y <- as.matrix(Y)
  tY <- as.matrix(tY)
  return(list(y = Y, ty = tY, train_index = train_index, test_index = test_index))
}

make_sparse_linear_model_cc <- function(n, p, h, s, rho, u, n0, h0, K, geno, model_name = "slm_cc") {
  data <- get_case_control_samples(n, geno, p, h, rho, s, u)
  beta <- data$beta
  N <- data$N
  se <- 4
  v <- 1 - (1 - u) * s / (1 - s)
  x <- data$x
  mu <- data$mu
  new_model(name = model_name, label = sprintf("n = %s, p = %s, h = %s, s = %s, rho = %s,  n0 = %s, u=%s, v=%s, K = %s, geno =%s", n, p, h, s, rho, n0, u, v, K,geno),
            params = list(x = x, mu = mu, beta = beta, n = n, N = N, cc = 1,
                          p = p, h = h, n0 = n0, se = se, geno = geno),
            simulate = function(x, se, mu, N, nsim) {
              out <- matrix(0, 1, nsim)
              out <- split(out, col(out))
              for(i in 1:nsim){
                PVal <- matrix(0, p, K)
                for(k in 1:K){
                  obj <- generate_pvalue(n0, p, h0, rho, (beta != 0), u, v,geno)
                  PVal[,k] <- as.matrix(obj$p)
                }
                simu_obj <- simulate_y(mu, se, n, N)
                obj <- SingleVariableScan(x[simu_obj$train_index,], simu_obj$y)
                draw <- new("Draw", y = simu_obj$y, pvalue = PVal, pvalue0 = as.matrix(obj$pval), ty = simu_obj$ty,
                            train_index = simu_obj$train_index, test_index = simu_obj$test_index)
                out[i] =  draw
              }
              return(out)
            })
}

make_sparse_linear_model <- function(n, p, h, s, rho, u, n0, h0, K, geno, model_name = "slm") {
  if(geno == 0){
    x <- createX(n, p, rho)
    tx <- createX(n, p, rho)
  }else{
    maf = 0.05 + 0.45*matrix(runif(p),1,p);
    x <- createGeno(n, p, rho, maf)
    tx <- createGeno(n, p, rho, maf)
  }
  v <- 1 - (1 - u) * s / (1 - s)
  se <- 4
  beta <- getbeta(x, p, h, s, se = se)
  mu <- as.numeric(x %*% beta)
  mut <- as.numeric(tx %*% beta)
  new_model(name = model_name, label = sprintf("n = %s, p = %s, h = %s, s = %s, rho = %s,  n0 = %s, u=%s, v=%s, K = %s, geno =%s", n, p, h, s, rho, n0, u, v, K,geno),
            params = list(x = x, beta = beta, mu = mu, se = se, n = n,
                          p = p, h = h, tx = tx, mut = mut, n0 = n0, geno = geno, maf = maf),
            simulate = function(x, mu, se,  mut, nsim) {
              out <- matrix(0, 1, nsim)
              out <- split(out, col(out))
              for(i in 1:nsim){
                PVal <- matrix(0, p, K)
                for(k in 1:K){
                  pvalue <- generate_pvalue(n0, p, h0, rho, (beta != 0), u, v,geno,maf)
                  PVal[,k] <- as.matrix(pvalue)
                }
                y <- as.matrix(mu + sqrt(se) * rnorm(n))
                ty <- as.matrix(mut + sqrt(se) * rnorm(n))
                obj <- SingleVariableScan(x, y)
                draw <- new("Draw", y = y, pvalue = PVal, pvalue0 = as.matrix(obj$pval), ty = ty)
                out[i] =  draw
              }
              return(out)
            })
}

make_sparse_linear_model_dependent <- function(n, p, h, s, rho, u, n0, h0, K, geno, uK, K0 = 0, model_name = "slm_dependent") {
  if(geno == 0){
    x <- createX(n, p, rho)
    tx <- createX(n, p, rho)
  }else{
    maf = 0.05 + 0.45*matrix(runif(p),1,p);
    x <- createGeno(n, p, rho, maf)
    tx <- createGeno(n, p, rho, maf)
  }
  v <- 1 - (1 - u) * s / (1 - s)
  vK <- 1 - (1 - uK) * s / (1 - s)
  se <- 4
  beta <- getbeta(x, p, h, s, se = se)
  mu <- as.numeric(x %*% beta)
  mut <- as.numeric(tx %*% beta)
  new_model(name = model_name, label = sprintf("n = %s, p = %s, h = %s, s = %s, rho = %s,  n0 = %s, u=%s, v=%s, K = %s, geno =%s, uK = %s", n, p, h, s, rho, n0, u, v, K,geno, uK),
            params = list(x = x, beta = beta, mu = mu, se = se, n = n,
                          p = p, h = h, tx = tx, mut = mut, n0 = n0, geno = geno, maf = maf, K0 = K0, K = K),
            simulate = function(x, mu, se,  mut, nsim) {
              out <- matrix(0, 1, nsim)
              out <- split(out, col(out))
              for(i in 1:nsim){
                PVal <- matrix(0, p, K + K0)
                obj <- generate_pvalue(n0, p, h0, rho, (beta != 0), u, v, geno, maf)
                PVal[,1] <- as.matrix(obj$p)
                gamma0 <- obj$gamma
                if(K >= 2){
                  for(k in 2:K){
                    # print(paste0("u=",u," v=",v," nz = ", sum(beta != 0)));
                    obj <- generate_pvalue(n0, p, h0, rho, gamma0, uK, vK, geno, maf)
                    PVal[,k] <- as.matrix(obj$p)
                  }
                }
                if(K0 >= 1){
                  for ( k in (1:K0) ) {
                    obj <- generate_pvalue(n0, p, h0, rho, gamma0, -1, -1, geno, maf)
                    PVal[,K + k] <- as.matrix(obj$p)
                  }
                }
                y <- as.matrix(mu + sqrt(se) * rnorm(n))
                ty <- as.matrix(mut + sqrt(se) * rnorm(n))
                obj <- SingleVariableScan(x, y)
                draw <- new("Draw", y = y, pvalue = PVal, pvalue0 = as.matrix(obj$pval), ty = ty, gamma_one = gamma0)
                out[i] =  draw
              }
              return(out)
            })
}

make_sparse_linear_model_dependentP_withSigma <- function(n, p, h, s, rho, u, v, n0, h0, sigma, geno , model_name = "slm_dependent_with_Sigma") {
  if(geno == 0){
    x <- createX(n, p, rho)
    tx <- createX(n, p, rho)
  }else{
    maf = 0.05 + 0.45*matrix(runif(p),1,p);
    x <- createGeno(n, p, rho, maf)
    tx <- createGeno(n, p, rho, maf)
  }
  se <- 4
  K = nrow(sigma)
  G <- simuGamma2(p, K + 1, dr, sigma)
  beta <- getbeta_give_gamma(x, p, h, G[,1], se)
  mu <- as.numeric(x %*% beta)
  mut <- as.numeric(tx %*% beta)
  new_model(name = model_name, label = sprintf("n = %s, p = %s, h = %s, s = %s, rho = %s,  n0 = %s, u=%s, v=%s, dr = %s", n, p, h, s, rho, n0, u, v, dr),
            params = list(x = x, beta = beta, mu = mu, se = se, n = n,
                          p = p, h = h, tx = tx, mut = mut, n0 = n0, G = G, maf = maf),
            simulate = function(x, mu, se,  mut, nsim, G = G) {
              out <- matrix(0, 1, nsim)
              out <- split(out, col(out))
              for(i in 1:nsim){
                PVal <- matrix(0, p, K)
                for(k in 1:K){
                  obj <- generate_pvalue(n0, p, h0, rho, (G[, k + 1] != 0), 1, 1, 1, maf)
                  PVal[,k] <- as.matrix(obj$p)
                }
                y <- as.matrix(mu + sqrt(se) * rnorm(n))
                ty <- as.matrix(mut + sqrt(se) * rnorm(n))
                obj <- SingleVariableScan(x, y)
                draw <- new("Draw", y = y, pvalue = PVal, pvalue0 = as.matrix(obj$pval), ty = ty,
                            train_index = 0, test_index = 0)
                out[i] =  draw
              }
              return(out)
            })
}



make_sparse_linear_model_dependent_pvalue <- function(n, p, h, s, rho, u, v, n0, h0, dr, geno , K = 1, model_name = "slm_dependent") {
  if(geno == 0){
    x <- createX(n, p, rho)
    tx <- createX(n, p, rho)
  }else{
    maf = 0.05 + 0.45*matrix(runif(p),1,p);
    x <- createGeno(n, p, rho, maf)
    tx <- createGeno(n, p, rho, maf)
  }
  se <- 4
  G <- simuGamma(p, K + 1, dr)
  beta <- getbeta_give_gamma(x, p, h, G[,1], se)
  mu <- as.numeric(x %*% beta)
  mut <- as.numeric(tx %*% beta)
  new_model(name = model_name, label = sprintf("n = %s, p = %s, h = %s, s = %s, rho = %s,  n0 = %s, u=%s, v=%s, dr = %s", n, p, h, s, rho, n0, u, v, dr),
            params = list(x = x, beta = beta, mu = mu, se = se, n = n,
                          p = p, h = h, tx = tx, mut = mut, n0 = n0, G = G, maf = maf),
            simulate = function(x, mu, se,  mut, nsim, G = G) {
              out <- matrix(0, 1, nsim)
              out <- split(out, col(out))
              for(i in 1:nsim){
                PVal <- matrix(0, p, K)
                for(k in 1:K){
                  obj <- generate_pvalue(n0, p, h0, rho, (G[, k + 1] != 0), 1, 1, 1, maf)
                  PVal[,k] <- as.matrix(obj$p)
                }
                y <- as.matrix(mu + sqrt(se) * rnorm(n))
                ty <- as.matrix(mut + sqrt(se) * rnorm(n))
                obj <- SingleVariableScan(x, y)
                draw <- new("Draw", y = y, pvalue = PVal, pvalue0 = as.matrix(obj$pval), ty = ty,
                            train_index = 0, test_index = 0)
                out[i] =  draw
              }
              return(out)
            })
}
