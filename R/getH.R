#' Get.H
#' @title get.H solve equation for H
#' @param W  weight matrix
#' @param y observations
#' @param mode This refers to the mode of resolution of matrix H. The default is "glmnet", user can chose "lars".
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#' @examples
#' p <- 5
#' n <- 10
#' W <- rbind(c(1,0,0,0,0),c(1,0,0,0,0),
#'           c(0,1,0,0,0),c(0,1,0,0,0) ,
#'            c(0,0,1,0,0),c(0,0,1,0,0) ,
#'             c(0,0,0,1,0),c(0,0,0,1,0) ,
#'              c(0,0,0,0,1),c(0,0,0,0,1) )
#' H <- as.matrix(rSpMatrix(p, 100, nnz = 100,
#' rand.x= function(nnz) round(rnorm(nnz, 5, 1), 2) ))
#' roc <- sapply(c(0.0,0.1, 0.3,0.5, 0.7, 1, 1.2,1.5), function(ss){
#'        y <- as.matrix(W%*%H + rnorm(20, sd=ss))
#'        test <- get.H(W, y=y, t="none", verbose=TRUE)
#'        tpr <- length(intersect(which(as.matrix(H)!=0),which(test!=0)))/sum(H!=0)
#'        fpr <- length(intersect(which(as.matrix(H)==0),which(test!=0)))/sum(H==0)
#'        return(c(fpr=fpr, tpr=tpr))
#'        })
#'
#'
#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#' @importFrom dplyr %>%
#' @export
get.H <- function(W, y, flavor_mod="glmnet", verbose) {
  J <- ncol(y)
  p <- ncol(W)
  if (J >1000){
    J1 <- 1000
  }else{J1 <- J}


  if(flavor_mod=="glmnet"){
    ## Lasso regression
    W.tilde <- Matrix::Matrix(kronecker(Matrix::Diagonal(J), W), sparse = TRUE)
    #W.tilde <- kronecker(diag(J), W) %>% as.big.matrix()
    y.tilde <- (as.numeric(y))
    cvfit <- try(cv.glmnet(W.tilde, y.tilde,
                           intercept=FALSE,
                           lambda = seq(from=0.1/J1, to =2/J1, length=20),
                           nfold=5))
    if (inherits(cvfit, "try-error")) {
      message("Something wrong occurs, init H to random")
      Z <-  matrix(0, nrow = p, ncol = J, byrow = FALSE)
      if(verbose) message(sprintf("number of non zero is %s on %s",sum(Z!=0), ncol(Z)*nrow(Z)))

    } else{
      Z <- as.matrix(coef(cvfit)[-1,]) %>% matrix(nrow = p, ncol = J, byrow = FALSE)
      if(verbose) message(sprintf("number of non zero is %s on %s",sum(Z!=0), ncol(Z)*nrow(Z)))
    }

    #system.time(cvfit <- cv.biglasso(W.tilde , y.tilde, seed = 1234,  nfolds=10, ncores=3, nlambda=10))
    #coefs <- as.matrix(coef(cvfit))
    ## Go back to Z
    #ind.best <- which(z.tilde$lambda==cvfit$lambda.1se)
    #Z <- as.matrix(z.tilde$beta)[, ind.best] %>% matrix(nrow = p, ncol = J, byrow = FALSE)
    res <- round(Z,2)
  }
  if(flavor_mod=="mglmnet"){
    cvfit <- try(cv.glmnet(W, y,family="mgaussian",
                           intercept=FALSE,
                           nfold=5, parallel=TRUE))
    Z <- as.matrix(do.call(cbind, coef(cvfit)) [-1,])


    #system.time(cvfit <- cv.biglasso(W.tilde , y.tilde, seed = 1234,  nfolds=10, ncores=3, nlambda=10))
    #coefs <- as.matrix(coef(cvfit))
    ## Go back to Z
    #ind.best <- which(z.tilde$lambda==cvfit$lambda.1se)
    #Z <- as.matrix(z.tilde$beta)[, ind.best] %>% matrix(nrow = p, ncol = J, byrow = FALSE)
    res <- round(Z,2)
  }
  print(flavor_mod)
  if(flavor_mod=="parallel"){
    if(verbose) message(sprintf("Parallele future Mode thanks HB"))

    Z <- future.apply::future_apply(y, 2,function(yy) {
      h <- cv.glmnet(x=W, y=yy, intercept = FALSE)
      coef(h)[-1,]
    })
    if(verbose) message(sprintf("number of non zero is %s on %s",sum(Z!=0), ncol(Z)*nrow(Z)))
    res <- round(Z,2)
  }
  return(res)
}
