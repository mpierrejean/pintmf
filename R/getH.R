#' Get.H
#' @title get.H solve equation for H
#' @param W  weight matrix
#' @param y observations
#' @param flavor_mod This refers to the mode of resolution of matrix H. The default is "glmnet", user can chose "oem" (not yet implemented).
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom ncvreg cv.ncvreg
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#' @importFrom future.apply future_apply
#' @importFrom dplyr %>%
#' @importFrom biglasso cv.biglasso
#' @importFrom methods slot
#' @importFrom quadrupen crossval
#' @importFrom bigmemory as.big.matrix
#' @export
get.H <- function(W, y, flavor_mod="glmnet", verbose) {
  J <- ncol(y)
  p <- ncol(W)
  if (J >1000){
    J1 <- 1000
  }else{J1 <- J}


  if(flavor_mod=="glmnet"){
    ## Lasso regression
    W.tilde <- Matrix(kronecker(Diagonal(J), W), sparse = TRUE)
    #W.tilde <- kronecker(diag(J), W) %>% as.big.matrix()
    y.tilde <- (as.numeric(y))
    cvfit <- try(cv.glmnet(W.tilde, y.tilde,
                           intercept=FALSE,
                           lambda = seq(from=0.1/J1, to =2/J1, length=20),
                           nfolds=5))
    if (inherits(cvfit, "try-error")) {
      message("Something wrong occurs, init H to random")
      Z <-  matrix(0, nrow = p, ncol = J, byrow = FALSE)
      if(verbose) message(sprintf("number of non zero is %s on %s",sum(Z!=0), ncol(Z)*nrow(Z)))

    } else{
      Z <- as.matrix(coef(cvfit)[-1,]) %>% matrix(nrow = p, ncol = J, byrow = FALSE)
      if(verbose) message(sprintf("number of non zero is %s on %s",sum(Z!=0), ncol(Z)*nrow(Z)))
    }
    Z <- round(Z,2)
  }
  if(flavor_mod=="ncvreg"){
    ## Lasso regression
    W.tilde <- Matrix(kronecker(Diagonal(J), W), sparse = FALSE) %>% as.matrix
    #W.tilde <- kronecker(diag(J), W) %>% as.big.matrix()
    y.tilde <- (as.numeric(y))
    cvfit <- try(cv.ncvreg(W.tilde, y.tilde,
                           intercept=FALSE,
                        penalty = "lasso",
                           nfolds=5))
    Z <-  coef(cvfit)[-1]%>% matrix(nrow = p, ncol = J, byrow = FALSE)
    Z <- round(Z,2)
  }
  if(flavor_mod=="biglasso"){
    ## Lasso regression
    W.tilde <- Matrix(kronecker(Diagonal(J), W), sparse = FALSE) %>% as.matrix%>% as.big.matrix()
    #W.tilde <- kronecker(diag(J), W) %>% as.big.matrix()
    y.tilde <- (as.numeric(y))
    fit <- cv.biglasso(X = W.tilde, y = y.tilde)
    Z <-  coef(fit)[-1]%>% matrix(nrow = p, ncol = J, byrow = FALSE)
  }
  if(flavor_mod=="quadrupen"){
    W.tilde <- Matrix(kronecker(Diagonal(J), W), sparse = FALSE) %>% as.matrix
    #W.tilde <- kronecker(diag(J), W) %>% as.big.matrix()
    y.tilde <- (as.numeric(y))
    beta.lasso <- slot(crossval(x=W.tilde,y=y.tilde, penalty="elastic.net", mc.cores=2,intercept=FALSE, K = 5) , "beta.min")
    Z <-   beta.lasso%>% matrix(nrow = p, ncol = J, byrow = FALSE)
  }
  return(Z)
}
