#' Get.H.group
#' @title get.H.group solve equation for H
#' @param W  weight matrix
#' @param y observations
#' @param type type of data (by default all are "none"),
#'   can take values into \code{ c("none", "methylation", "mutation")},
#'   this argument allows to adapt modelisation by adding constraints.
#'   For exemple methylation data takes values between 0 and 1, then H block corresponding to
#'   methylation data takes values between 0 and 1.
#' @param group a vector of the same length of y
#' @examples
#' library(tidyverse)
#' p <- 5
#' n <- 10
#' W <- buildWeightMatrix(n,p)
#' n.vars <- 100
#' H <- rbind(c(runif(15,  0, 5), rep(0, n.vars - 15)),
#'            c(rep(0, n.vars - 25), runif(10, 0, 5), rep(0, n.vars - 85)),
#'            c(rep(0, n.vars - 50), runif(20,  -5, 0), rep(0, n.vars - 70)),
#'            c(rep(0, n.vars - 40), runif(15,  -5, 0), rep(0, n.vars - 75)),
#'            c(rep(0, n.vars - 70), runif(20, 0, 5), rep(0, n.vars - 50)))  %>%
#'            round(2)
#' groups = c(rep(1, 15), rep(2, 15), rep(3, 20), rep(4, 10), rep(5, 10), rep(6,5), rep(7,10), rep(8, 15))
#' y <- as.matrix(W%*%H + rnorm(20, sd=0.0))
#' test <- get.H.group(W, y=y, t="none", groups=groups)
#' tpr <- length(intersect(which(as.matrix(H)!=0),which(test!=0)))/sum(H!=0)
#' tpr
#' fpr <- length(intersect(which(as.matrix(H)==0),which(test!=0)))/sum(H==0)
#' fpr
#' @importFrom oem oem
#' @importFrom dplyr %>%
get.H.group <- function(W, y, t, groups) {
  #stopifnot(length(lambda) == 1L)  ## sanity check
  J <- ncol(y)
  p <- ncol(W)
  ## Lasso regression
  W.tilde <- Matrix::Matrix(kronecker(Matrix::Diagonal(J), W), sparse = TRUE)
  y.tilde <- (as.numeric(y))

  ##  z.tilde <- scoop::lasso(W.tilde,y.tilde ,intercept=FALSE, n.lambda=20,
  ## lambda.min=lambda, verbose=TRUE, normalize=FALSE)
  ##  z.tilde <- lars::lars(W.tilde,y.tilde ,type="lar", intercept=FALSE,
  ##      normalize=FALSE)
  if (t == "mutation" || t == "methylation" ) {
    ll  <- 0
  }else{
    ll = -Inf
  }
  if (t == "methylation") {
    ul  <- 1
  }else{
    ul = Inf
  }
  grp.lam <- oem(W.tilde, y.tilde, penalty = "sparse.grp.lasso",
                 groups = rep(groups,p), intercept=FALSE)
  idx <- which(grp.lam$nzero[[1]]!=0)
  loss <- apply(grp.lam$beta$sparse.grp.lasso[,idx],2, function(beta){
    beta2 <- beta[-1]
    Z <- as.matrix(beta2) %>% matrix(nrow = p, ncol = J, byrow = FALSE)
    sum((y-W%*%Z)^2)/(J*nrow(y))
  })
  best.mod <- which(abs(diff(loss))<1e-3) %>% min
  ## Go back to Z
  beta <- grp.lam$beta$sparse.grp.lasso[-1,best.mod]
  Z <- as.matrix(beta) %>% matrix(nrow = p, ncol = J, byrow = FALSE)
   return(round(Z,2))
}
