context("Consistency of the optimization algorithm: get.W")
set.seed(33)
W <- rbind(diag(1, nrow = 3, ncol = 3), diag(1, nrow = 3, ncol = 3))
H1 <- as.matrix(rSpMatrix(3, 100, nnz = 10, rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) ))
H2 <- as.matrix(rSpMatrix(3, 100, nnz = 10,  rand.x= function(nnz) round(rnorm(nnz, 2, 0.2), 2) ))
H3 <- as.matrix(rSpMatrix(3, 300, nnz = 10, rand.x= function(nnz) round(rnorm(nnz, 6, 0.2), 2) ))
Y1 <- as.matrix(W%*%H1 + rnorm(20, sd=1))
Y2 <- as.matrix(W%*%H2 + rnorm(20, sd=1))
Y3 <- as.matrix(W%*%H3 + rnorm(20, sd=1))
Ybar <- cbind(Y1, Y2, Y3)
Zbar <- cbind(H1, H2, H3) %>% t
group <- apply(W, 1, which.max)
test_that("Consistency of get.W", {
  ## solving with PintMF
  We <- PintMF::get.W(Zbar, Ybar)
  print(We)
  expect_equal(is.matrix(We), TRUE)
  expect_equal(ncol(We), ncol(W))
  expect_equal(nrow(We), nrow(W))
})



test_that("Consistency of get.W in supervised way", {
  ## solving with PintMF
  We <- PintMF::get.W.supervised(Zbar, Ybar, group=group)
  expect_equal(is.matrix(We), TRUE)
  expect_equal(ncol(We), ncol(W))
  expect_equal(nrow(We), nrow(W))
  print(We)
})
