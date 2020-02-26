context("Consistency of the optimization algorithm: get.H")
set.seed(33)
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



