library(simulator)

fdr2FDR <- function(fdr){
  M <- length(fdr)
  sortfdr <- sort.int(fdr, index.return=TRUE)
  sfdr <- sortfdr$x
  index <- sortfdr$ix
  gFDR <- cumsum(sfdr) / as.matrix(1:M, M,1)
  gFDR[gFDR > 1] = 1;
  gFDR[index] = gFDR;
  return(gFDR)
}

sqrerr <- new_metric("MSE", "Mean Squared Error for Beta Estimation",
                     metric = function(model, draw, out) {
                       if( "beta" %in% names(out)){
                         colMeans(as.matrix(out$beta - model$beta)^2)
                       }else{
                          NA
                        }

                     })

corr <- new_metric("CORR", "Prediction Accuracy",
                  metric = function(model, draw, out) {
                    if(sum(out$tyhat) != 0 )
                      result <- cor(out$tyhat, draw@ty)
                    else
                      result <- -1
                  })

power <- new_metric("Power", "Power of Risk Variants Identification",
                   metric = function(model, draw, out) {
                    if( "fdr" %in% names(out))
                    {
                       gFDR <- fdr2FDR(out$fdr)
                       found <- gFDR < 0.1
                       nzeros <- (model$beta != 0)
                       sum(nzeros[found] == 1) / sum(nzeros)
                    }else{
                      NA
                    }
                   })

rauc <- new_metric("AUC", "Risk Variant Identification ",
                   metric = function(model, draw, out) {

                     if( "fdr" %in% names(out))
                    {
                       nzeros <- (model$beta != 0)
                       pROC::auc(nzeros, out$fdr)
                     }else{
                       NA
                     }

                   })


cauc <- new_metric("CAUC", "Classification Accuracy",
                   metric = function(model, draw, out) {
                     if( "tyhat" %in% names(out))
                     {
                       pROC::auc(as.numeric(draw@ty),as.numeric(out$tyhat))
                     }else{
                       NA
                     }
                   })

FDR <- new_metric("FDR", "False Discovery Rate",
                  metric = function(model, draw, out) {
                    result <- NA
                    if( "fdr" %in% names(out))
                    {
                      FDR <- fdr2FDR(out$fdr)
                      nzeros <- (model$beta != 0)
                      found <- FDR < 0.1
                      nfound <- sum(found)
                      errVE <- sum(nzeros[found] == 0)
                      if(nfound != 0){
                        result = errVE / nfound
                      }else{
                        result = 0
                      }
                    }
                    result
                  })

R2 <- new_metric("R2", "prediction accuracy",
                 metric = function(model, draw, out) {
                   # cor(out$tyhat, draw@ty)
                   div = draw@ty - mean(draw@ty)
                   diff =  draw@ty - out$tyhat;
                   1 - sum(diff * diff) / (sum(div * div));
                 })


sqrerr <- new_metric("sqrerr", "beta squared error",
                     metric = function(model, draw, out) {
                       if( "beta" %in% names(out)){
                         colMeans(as.matrix(out$beta - model$beta)^2)
                       }else{
                         NA
                       }

                     })


mse <- new_metric("mse", "y mean square error",
                  metric = function(model, draw, out) {
                    # result <- "NA"
                    # if(!is.null(out$tyhat))
                    # {
                    if(sum(out$tyhat) != 0 )
                      result <- sqrt(colMeans(as.matrix((out$tyhat - draw@ty))^2))
                    else
                      result <- -1
                    # }
                    # result

                  })


mae <- new_metric("mae", "mean absolute error",
                   metric = function(model, draw, out) {
                     # result <- "NA"
                     # if(!is.null(out$tyhat))
                     # {
                     if(sum(out$tyhat) != 0 )
                       result <- colMeans(as.matrix(abs(out$tyhat - draw@ty)))
                     else
                       result <- -1
                     # }
                     # result

                   })


