#' computeLL
#' @title computeLL Compute likelihood of the model
#' @param Y observations
#' @param res result from SolveInt function
#' @return a list containing RSS value, p, n and J
#' @export
computeLL <- function(Y, res){
  K <- length(Y)
  W <- res$W
  p <- ncol(W)
  n <- nrow(Y[[1]])
  sub_rss_k <- sapply(1:K, function(kk){
    yk <- Y[[kk]]
    Jk <- ncol(yk)
    Hk <- res$H[[kk]]
    sum((yk-W%*%Hk)^2)/(n*Jk)
  })
  J <- sum(sapply(Y, ncol))
  RSS <- sum(sub_rss_k)/K
  return(c(RSS=RSS, p=p, n=n, J=J))
}

#' compute_BIC
#' @title compute_BIC Compute BIC the model
#' @param LL result from computeLL function
#' @param pen penalty value in formula
#' @export
compute_BIC <- function(LL, pen=1){
  BIC <-LL["n"]*log(LL["RSS"])+pen*LL["p"]*log(LL["n"])
  return(c(BIC=BIC, p=LL["p"] ))
}

#' compute_coph
#' @title compute_coph the cophenetic correlation coefficient
#' @param res result from SolveInt function
#'
#' @return coph the cophenetic correlation coefficient
#' @export
compute_coph <- function(res){
  coph <- cor(hclust(dist(res$W), method="ward.D2") %>% cophenetic, dist(res$W))
  return(coph)
}
