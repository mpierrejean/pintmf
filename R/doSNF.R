#' init_SNF
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param sigma Variance for local model
#' @param K_n Number of nearest neighbors
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method SNF
#' @import SNFtool
#' @importFrom dplyr %>%
#' @export
init_SNF <- function (data, K,  K_n=10, sigma=0.5){
  dat <- lapply(data, function (dd){
    dd <- dd %>% as.matrix
    W <- dd %>% dist2(dd) %>% affinityMatrix(K=K_n, sigma=sigma)
  })
  W <-  SNF(dat, K_n, K_n)
  clust.SNF = W %>% spectralClustering(K)
  res <- list(clust=clust.SNF, fit= W)
  return(res)
}
