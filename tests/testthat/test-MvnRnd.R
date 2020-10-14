library(MomentMatching)

rets  <- diff(log(EuStockMarkets))
mu    <- colMeans(rets)
sigma <- stats::cov(rets)
simul <- MvnRnd(M = mu, S = sigma, J = 50)

test_that("J argument rejects vectors", {
  expect_error(MvnRnd(M = mu, S = sigma, J = c(1, 1)))
})

test_that("simul has the right dimensions", {
  expect_equal(nrow(simul), 50)
  expect_equal(ncol(simul), 4)
})

test_that("the output is a matrix", {
  expect_true(is.matrix(simul))
})

test_that("Transposes the argument 'mu' when a nx1 matrix is provided", {
  expect_equal(nrow(MvnRnd(M = t(t(mu)), S = sigma, J = 50)), 50)
})
