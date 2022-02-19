context("poisthin")

test_that("poisthin input works", {
  ## Generate example data ---------------------------------------------------
  n <- 100
  p <- 1000
  nsamp <- 10
  ngene <- 100
  mat <- matrix(stats::rpois(n = n * p, lambda = 20), nrow = n)
  signal_fun <- stats::rnorm
  signal_params <- list(mean = 0, sd = 1)
  gselect <- "max"
  gvec <- rep(FALSE, length = ncol(mat))
  gvec[sample(1:ncol(mat), size = ngene)] <- TRUE
  prop_null <- 0.9

  args <- list(mat = mat, nsamp = nsamp, ngene = ngene,
               gselect = "max",
               gvec = NULL,
               skip_gene = 1,
               signal_fun = stats::rnorm,
               signal_params = list(mean = 0, sd = 1),
               prop_null = 0.5)

  pout <- do.call(what = poisthin, args = args)
  expect_equal(mean(pout$beta == 0), 0.5)
  expect_equal(ngene, ncol(pout$Y))
  expect_equal(nsamp, nrow(pout$Y))

  args$gselect <- "rand_max"
  pout <- do.call(what = poisthin, args = args)
  expect_equal(ngene, ncol(pout$Y))
  expect_equal(nsamp, nrow(pout$Y))

  args$gselect <- "random"
  pout <- do.call(what = poisthin, args = args)
  expect_equal(ngene, ncol(pout$Y))
  expect_equal(nsamp, nrow(pout$Y))

}
)

test_that("poisthin returns about the correct effect sizes", {
  set.seed(15)
  mat <- matrix(stats::rpois(n = 10000, lambda = 1000), ncol = 1)
  true_esize <- 2
  signal_fun <- function(n, true_esize) { rep(true_esize, n) }
  pout <- poisthin(mat = mat, nsamp = nrow(mat), ngene = 1,
                   signal_fun = signal_fun,
                   signal_params = list(true_esize = true_esize),
                   prop_null = 0 )


  ly <- log(pout$Y + 1, base = 2)
  esize <- mean(ly[pout$X[, 2] == 0], na.rm = TRUE) - mean(ly[pout$X[, 2] == 1], na.rm = TRUE)

  expect_equal(esize, true_esize, tolerance = 0.002)
}
)



test_that("poisthin input works with many zeros", {
  ## Generate example data ---------------------------------------------------
  n <- 100
  p <- 1000
  nsamp <- 10
  ngene <- 100
  mat <- matrix(stats::rpois(n = n * p, lambda = 1), nrow = n)
  signal_fun <- stats::rnorm
  signal_params <- list(mean = 0, sd = 1)
  gselect <- "max"
  gvec <- rep(FALSE, length = ncol(mat))
  gvec[sample(1:ncol(mat), size = ngene)] <- TRUE
  prop_null <- 0.9

  args <- list(mat = mat, nsamp = nsamp, ngene = ngene,
               gselect = "max",
               gvec = NULL,
               skip_gene = 1,
               signal_fun = stats::rnorm,
               signal_params = list(mean = 0, sd = 1),
               prop_null = 1)

  expect_error(pout <- do.call(what = poisthin, args = args), NA)

}
)
