#' PVE computing
#'
#' @param data alist of original data of size $n x J_k$
#' @param H a list of feature matrix $p x J_k$
#' @param W a matrix of weights size $n x p$
#'
#' @return pve percentage of variance explained
#' @export
#'
PVE <- function(data, H, W){
  ind <- 1:length(data)
  pve <- 100*mean(sapply(ind, function(ss){
    Y <- data[[ss]]
    h <- H[[ss]]
    Ybar <- rowMeans(Y)
    pve_s <- 1-sum((Y-W%*%h)^2)/sum((Y-Ybar)^2)
    return(pve_s)
  }))
  return(pve)
}
