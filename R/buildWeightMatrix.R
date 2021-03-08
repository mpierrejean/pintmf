#' Create a simulated weight matrix
#'
#' @param nrow An integer value, the number of patients (in row)
#'
#' @param ncol An integer value, the number of clusters (in columns)
#'
#' @param prop An numerical value between 0 and 1, the proportion of non-zero coefficient
#' by row, by default 0.5
#'
#' @return A sparse \code{n} x \code{p} matrix
#'
#' @examples
#' W <- buildWeightMatrix(10,6)
#' @importFrom Matrix spMatrix
#' @export
buildWeightMatrix <- function(nrow, ncol, prop=0.5) {
  nnz <- round(prop*ncol)
  rand.x = function(n) runif(n,0,5)
  W <- matrix(0,nrow, ncol)
  while(!all(colSums(W)>0)|!all(rowSums(W)>0)  ## sanity check
){
    cat("starting simulation\n")
    W <- do.call(rbind, replicate(nrow,
                                  spMatrix(1, ncol,
                                                         i = sample(1, nnz, replace = TRUE),
                                                         j = sample(ncol, nnz, replace = TRUE),
                                                         x = rand.x(nnz))))
    W <- as.matrix(W)
  }
  cat("done")
  W <- round(W/rowSums(W),2)
  return(W)
}
