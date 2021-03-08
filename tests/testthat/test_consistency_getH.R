context("Consistency of the optimization algorithm: get.H")
set.seed(33)

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
  spMatrix(nrow, ncol,
           i = sample(nrow, nnz, replace = TRUE),
           j = sample(ncol, nnz, replace = TRUE),
           x = rand.x(nnz))
}
W <- rbind(diag(1, nrow = 3, ncol = 3), diag(1, nrow = 3, ncol = 3))
H1 <- as.matrix(rSpMatrix(3, 100, nnz = 10, rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) ))
Y1 <- as.matrix(W%*%H1 + rnorm(20, sd=1))
test_that("Consistency of get.W", {
  ## solving with PintMF
  He <- PintMF::get.H(W, Y1, flavor_mod = "glmnet", verbose=TRUE)
  expect_equal(is.matrix(He), TRUE)
  expect_equal(ncol(He), ncol(H1))
  expect_equal(nrow(He), nrow(H1))
  He <- PintMF::get.H(W, Y1, flavor_mod = "ncvreg", verbose=TRUE)
  expect_equal(is.matrix(He), TRUE)
  expect_equal(ncol(He), ncol(H1))
  expect_equal(nrow(He), nrow(H1))

})



