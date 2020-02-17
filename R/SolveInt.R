#' SolveInt
#' @title SolveInt solve equation for H and W
#'
#' @param Y observations, list of blocks of data
#' @param p The number of latent profiles
#' @param max.it by default 20. Maximum iteration of the algorithm, else until convergence
#'   of matrix W
#' @param type type of data (by default all are "none"),
#'   can take values into \code{ c("none", "methylation", "mutation")},
#'   this argument allows to adapt modelisation by adding constraints.
#'   For exemple methylation data takes values between 0 and 1, then H block corresponding to
#'   methylation data takes values between 0 and 1.
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#' @param flavor_mod glmnet or (oem not implemented yet)
#' @param init_flavor flavor to initialize W matrix (hclust, SNF, random, pca or svd)
#' @param ... other arguments
#'
#' @return A list that contains :
#'   \code{H} is a list of matrix (latent profiles for each block of data)
#'   \code{W} is a matrix (weight matrix).
#'   \code{loss} vector of the loss of W between two iterations
#'   \code{pve} Percentage of variance explained
#' @export
#' @importFrom glmnet glmnet
#' @importFrom dplyr bind_cols
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom nnet class.ind
#' @importFrom pcaMethods pca
#' @import stats
SolveInt <- function(Y, p, max.it=20, flavor_mod="glmnet",  type=rep("none", length(Y)), init_flavor="hclust", verbose=FALSE,...
) {
  if (!is.list(Y)) {
    stop("Y is not a list")
  }
  family = c("mutation", "none", "methylation")
  if (!all(type %in% family)) {
    stop(sprintf("Type is not in the list, choose type in %s", paste(family, collapse = ", ")))
  }

  if(verbose) message(sprintf("Run algorithm for p=%s\n", p))

  if(verbose) message(sprintf("type of initialization : %s\n", init_flavor))
  gaussian_type <- which(type == "none")
  Yt <- lapply(1:length(Y),
               function(yy) {
                 if (yy %in% gaussian_type) {
                   if(verbose) message (sprintf("Center cols of Y\n"))
                   y <- sweep(Y[[yy]],MARGIN = 2,colMeans(Y[[yy]]) )
                   # y <- scale(Y[[yy]])

                 }else{
                   y <- Y[[yy]]
                 }
                 return(y)
               })
  Yt <- Y
  dd <- lapply(Yt, dist)
  hc <- lapply(dd, hclust, method = "ward.D")
  ## initialization of H

  initH_clust <- function(yyt, hhc, pp) {
    cluster <- cutree(hhc, k = pp)
    ## Averaged profiles by cluster
    t(sapply(split(as.data.frame(yyt), f = cluster), FUN = colMeans))
  }
  initH_random <- function(yyt, pp) {
    idx <- sample(1:nrow(yyt),pp)
    yyt[idx,]
  }
  initH_pca <- function(yyt, p) {
    lapply(yyt, function (h) pca(h, nPcs = p)@loadings %>% t)
  }
  initH_svd <- function(yyt, p) {
    lapply(yyt, function (h) svd(h, nv=p)$v %>% t)
  }
  initH_snf <- function(yyt, p) {
    cluster <- init_SNF(yyt, K=p)$clust
    ## Averaged profiles by cluster
    lapply(yyt, function (y) t(sapply(split(as.data.frame(y), f = cluster), FUN = colMeans)))
  }
  if(verbose) message (sprintf("Init Hc \n"))

  if(init_flavor=="hclust"){
    Hc <- mapply(initH_clust,Yt,hc, rep(p, length(Y)), SIMPLIFY = FALSE)
  }
  if(init_flavor=="svd"){
    Hc <- initH_svd(Yt, p )
  }
  if(init_flavor=="SNF"){
    Hc <- initH_snf(Yt, p )
  }
  if(init_flavor=="pca"){
    Hc <- initH_pca(Yt, p )
  }else{
    Hc <- mapply(initH_random,Yt, rep(p, length(Y)), SIMPLIFY = FALSE)
  }

  it <- 1
  loss <- c(100)
  pve <- numeric(0)
  print(p)
  Wc <- matrix(0, nrow(Y[[1]]), p)
  while (it <= max.it) {
    if(verbose) message("Solve W\n")
    W.old <- Wc
    Hc_norm <- lapply(Hc, function(hh) hh/sqrt(ncol(hh)))
    Zbar <- t(bind_cols(lapply(Hc_norm, data.frame))) %>% as.matrix()
    Ybar <- bind_cols(lapply(Yt, data.frame)) %>% as.matrix
    if(it==1 & init_flavor=="snf"){
      Wc <- class.ind(init_SNF(Yt, K=p)$clust)
    }else{
      Wc <- get.W(Zbar=Zbar,
                  Ybar=Ybar)
    }
    if(verbose) message ("Solve Hc\n")
    Hc <- future_lapply(1:length(Yt), function(yy){
      if(verbose){
        message(sprintf("Solve data number %s:",yy))
      }
      print(flavor_mod)
      get.H(y = Yt[[yy]], W = Wc, flavor_mod=flavor_mod, verbose)}
    )


    loss[it + 1] <- sum((Wc - W.old)^2)/(ncol(Wc)*nrow(Wc))
    pve[it] <- PVE(Yt, Hc, Wc)
    if(verbose) (sprintf("loss : %s\n", round(loss[it + 1],2)))
    it <- it + 1
    if(verbose) message(sprintf("iteration number %s", cat(it, "\n")))
  }
  stop_ind=FALSE
  if (!all(colSums(Wc) > 0)) {
    cat("At least one column is zero, you should be run the algorithm with less latent variables\n")
    stop_ind=TRUE
  }
  ## Define names of variables
  var_names <- lapply(Y, colnames)
  Hc <- lapply(1:length(Y), function(ii) {
    h <- Hc[[ii]]
    colnames(h) = var_names[[ii]]
    return(h)
  })
  return(list(H = Hc, W = Wc, loss = loss, pve = pve, stop_ind=stop_ind))
}
