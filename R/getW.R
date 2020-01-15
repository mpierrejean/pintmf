#' Get.W
#' @title get.W solve equation of W
#' @param Zbar summary of various blocks of ceofficient  matrices
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
#' @importFrom lsei lsei
#' @importFrom future.apply future_lapply
#' @importFrom dplyr %>%
#' @export
get.W <- function(Zbar, Ybar) {
  plan(multiprocess)
  Wc <-  future_lapply(1:nrow(Ybar), FUN=solvew, Ybar=Ybar, Zbar=Zbar) %>% simplify2array %>% t
  return(Wc)
}


solvew <- function(ind, Ybar, Zbar){
  ybar <- Ybar[ind,]
  p <- ncol(Zbar)
  fit <-  try(glmnet(x=Zbar, y=ybar,
                      alpha=0,
                      lower.limit=0,lambda = 1,
                      intercept=FALSE))
  if (inherits(fit, "try-error")) {
    message("Something wrong occurs, init W with Lsei")
    w <- try(lsei::lsei(a=Zbar, b=ybar, c=rep(1, p), d=1,lower=0))
    if (inherits(w, "try-error")) {
      w <- rep(1/p, p)
    }
  } else{
     if(is.na(sum(as.vector(fit$beta)))){
       w <- rep(1/p, p)
     }else{
       w <- as.vector(fit$beta)
       #w <- w/sqrt(sum(w^2))
       w <- w/sum(w)
   # w <- fit
    }

  }
  return(round(w, 1))
}
