#' Get.W
#' @title get.W solve equation of W
#' @param Zbar summary of various blocks of ceofficient  matrices
#' @param flavor_mod This refers to the mode of resolution of matrix W. The default is "glmnet" else you can choose sparse_lsei,sparse_glmnet or cv_glmnet
#' @param Ybar summary of various block of  observations
#' library(Matrix)
#' library(dplyr)
#' W <- matrix(round(runif(20*2, 0, 5), 2), nrow=20, ncol=2)
#' W <- W/rowSums(W)
#' H1 <- as.matrix(rSpMatrix(2, 100000, nnz = 10000, rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) ))
#' H2 <- as.matrix(rSpMatrix(2, 1000, nnz = 100,  rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) ))
#' H3 <- as.matrix(rSpMatrix(2, 30000, nnz = 1000, rand.x= function(nnz) round(rnorm(nnz, 6, 0.2), 2) ))
#' Y1 <- as.matrix(W%*%H1 + rnorm(20, sd=1))
#' Y2 <- as.matrix(W%*%H2 + rnorm(20, sd=1))
#' Y3 <- as.matrix(W%*%H3 + rnorm(20, sd=1))
#' Ybar <- cbind(Y1, Y2, Y3)
#' Zbar <- cbind(H1, H2, H3) %>% t
#' We <- get.W(Zbar, Ybar)
#'
#' @importFrom glmnet glmnet
#' @importFrom limSolve lsei
#' @importFrom dplyr %>%
#' @export
get.W <- function(Zbar, Ybar, flavor_mod='glmnet') {
  Wc <-  lapply(1:nrow(Ybar), FUN=solvew, Ybar=Ybar, Zbar=Zbar, flavor_mod= flavor_mod) %>% simplify2array %>% t
  return(Wc)
}


solvew <- function(ind, Ybar, Zbar, flavor_mod){
  ybar <- Ybar[ind,]
  p <- ncol(Zbar)
  if(flavor_mod=='glmnet'){
    fit <-  try(glmnet(x=Zbar, y=ybar,
                       alpha=0,
                       lower.limits=0,lambda = 1,
                       intercept=FALSE))
    if (inherits(fit, "try-error")|sum(fit$beta)==0) {
      message("Something wrong occurs, init W with Lsei")
      w <- lsei(A=Zbar, B=ybar, H = rep(0, p), G = diag(1, p), type=2)$X
      w <- w/sum(w)

    } else{
      if(is.na(sum(as.vector(fit$beta)))){
        cvfit <- cv.glmnet(x=Zbar, y=ybar,
                           alpha=0,
                           lower.limits=0,
                           intercept=FALSE)
        w <- coef(cvfit)[-1,]
        w <- w/sum(w)
      }else{
        w <- as.vector(fit$beta)
        w <- w/sum(w)
      }

    }
  }
  else if(flavor_mod=='sparse_lsei'){
    w <- lsei(A=Zbar, B=ybar, H = rep(0, p), G = diag(1, p), type=2)$X
    w <- w/sum(w)
  }
  else if(flavor_mod=='sparse_glmnet'){
    w <- glmnet(x=Zbar, y=ybar,
                alpha=0,
                lower.limits=0,lambda = 0,
                intercept=FALSE)$beta %>% as.vector()
    w <- w/sum(w)
  }
  else if(flavor_mod=='cv_glmnet'){
    cvfit <- cv.glmnet(x=Zbar, y=ybar,
                       alpha=0,
                       lower.limits=0,
                       intercept=FALSE)
    w <- coef(cvfit)[-1,]
    w <- w/sum(w)
  }else{
    stop("Please choose flavor_mod : 'glmnet' or cv_glmnet or 'sparse_lsei','sparse_glmnet' or for matrix W")

  }
  return(round(w, 1))
}

