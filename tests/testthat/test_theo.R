context("Theoretical Model")


test_that("seqgen_base works", {
  set.seed(1)
  nsamp <- 100
  ngene <- 101

  designmat <- model.matrix(~condition,
                            data.frame(condition = factor(rep(c("A", "B"),
                                                              length.out = nsamp))))

  coefmat <- matrix(1, nrow = ngene, ncol = 2)
  dispvec <- rep(1, length.out = ngene)

  temp <- seqgen_base(designmat = designmat, coefmat = coefmat, dispvec = dispvec)
  meanvec <- colMeans(temp)
  expect_equal(mean(meanvec[designmat[, "conditionB"] == 1]), 4, tol = 0.1)
  expect_equal(mean(meanvec[designmat[, "conditionB"] == 0]), 2, tol = 0.1)

  designmat <- matrix(1, nrow = nsamp)
  coefmat   <- matrix(10, nrow = ngene)
  dispvec <- c(1, rep(0.00001, length.out = ngene - 1))
  temp <- seqgen_base(designmat = designmat, coefmat = coefmat, dispvec = dispvec)
  varvec  <- apply(temp, 1, var)
  expect_true(varvec[1] > max(varvec[-1]))

})


test_that("seqgen_diff works", {
  nsamp <- 10
  ngene <- 1000
  nvar  <- 2
  designmat <- matrix(rnorm(nvar * nsamp), nrow = nsamp)
  coefmat   <- matrix(rnorm(nvar * ngene), nrow = ngene)
  dispvec <- rchisq(n = ngene, df = 1)

  thout <- seqgen_diff(designmat = designmat,
                       coefmat   = coefmat,
                       dispvec   = dispvec)
  expect_equal(class(thout), "list")
})

test_that("seqgen_2group works", {
  ngene            <- 100
  nsamp            <- 10
  dispvec          <- rep(0.1, times = ngene)
  prop_null        <- 1
  signal_fun       <- stats::rnorm
  signal_params    <- list(mean = 0, sd = 1)
  intercept_fun    <- stats::rnorm
  intercept_params <- list(mean = 4, sd = 2)
  libvec           <- rep(1, times = nsamp)
  group_prop       <- 0.5
  design_sv        <- NULL
  coef_sv          <- NULL

  thout <- seqgen_2group(ngene = ngene, nsamp = nsamp)
  expect_equal(class(thout), "list")
})
