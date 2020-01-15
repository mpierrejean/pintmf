#' computeLL
#' @title computeLL Compute likelihood of the model
#' @param Y observations
#' @param res result from SolveInt function
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
    sum((yk-W%*%Hk)^2)/(Jk*p)
  })
  J <- sum(sapply(Y, ncol))
  RSS <- sum(sub_rss_k)
  return(c(RSS=RSS, p=p, n=n, J=J))
}

#' compute_BIC
#' @title compute_BIC Compute BIC the model
#' @param LL result from computeLL function
#' @param pen penalty value in formula
#' @export
compute_BIC <- function(LL, pen=1){
  BIC <- LL["n"]*log(LL["RSS"]/LL["n"])+pen*LL["p"]*log(LL["n"]*LL["J"])
  return(c(BIC=BIC, p=LL["p"] ))
}
