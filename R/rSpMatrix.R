#' rSpMatrix random sparse matrix
#' @title rSpMatrix create a non negative matrix
#' @param nrow number of row
#' @param ncol number of column
#' @param nnz sparsity coefficient
#' @param rand.x function to simulate data
#' @examples
#' library(Matrix)
#' H <- rSpMatrix(2, 500, nnz = 500, rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) )
#' @export
rSpMatrix <- function(nrow, ncol, nnz,
                      rand.x = function(nnz) round(rnorm(nnz), 2))
{
  ## Purpose: random sparse matrix
  ## --------------------------------------------------------------
  ## Arguments: (nrow,ncol): dimension
  ##          nnz  :  number of non-zero entries
  ##         rand.x:  random number generator for 'x' slot
  ## --------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14.-16. May 2007
  stopifnot((nnz <- as.integer(nnz)) >= 0,
            nrow >= 0, ncol >= 0, nnz <= nrow * ncol)
  Matrix::spMatrix(nrow, ncol,
                   i = sample(nrow, nnz, replace = TRUE),
                   j = sample(ncol, nnz, replace = TRUE),
                   x = rand.x(nnz))
}
