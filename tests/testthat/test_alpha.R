context("Alpha")

test_that("alpha = 1 is approximately correct", {
  set.seed(15)
  mat <- matrix(stats::rpois(n = 10000, lambda = 1000), ncol = 1)
  true_esize <- 2
  signal_fun <- function(n, true_esize) { rep(true_esize, n) }
  pout <- poisthin(mat = mat, nsamp = nrow(mat), ngene = 1,
                   signal_fun = signal_fun,
                   signal_params = list(true_esize = true_esize),
                   prop_null = 0, alpha = 1)

  sdhat <- apply(log2(mat + 1), 2, sd) / sqrt(nrow(mat))

  ly <- log2(pout$Y + 1)
  esize <- mean(ly[pout$X[, 2] == 0], na.rm = TRUE) - mean(ly[pout$X[, 2] == 1], na.rm = TRUE)
  expect_equal(esize, true_esize * sdhat, tol = 0.001)

})

test_that("alpha = 0.5 is approximately correct", {
  set.seed(15)
  mat <- matrix(stats::rpois(n = 10000, lambda = 1000), ncol = 1)
  true_esize <- 2
  signal_fun <- function(n, true_esize) { rep(true_esize, n) }
  pout <- poisthin(mat = mat, nsamp = nrow(mat), ngene = 1,
                   signal_fun = signal_fun,
                   signal_params = list(true_esize = true_esize),
                   prop_null = 0, alpha = 0.5)

  sdhat <- apply(log2(mat + 1), 2, sd) / sqrt(nrow(mat)) ## sd of the mean

  ly <- log2(pout$Y + 1)
  esize <- mean(ly[pout$X[, 2] == 0], na.rm = TRUE) - mean(ly[pout$X[, 2] == 1], na.rm = TRUE)
  expect_equal(esize, true_esize * sqrt(sdhat), tol = 0.001)

})
