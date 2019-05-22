context("general thinning")


test_that("general thinning works", {
  n <- 10
  p <- 100
  Z <- rnorm(n)
  alpha <- rnorm(p)
  mat <- round(2^(alpha %*% t(Z) + matrix(rnorm(n * p), nrow = p, ncol = n)))
  design_perm <- matrix(runif(n))
  coef_perm <- matrix(1, nrow = p, ncol = 1)
  target_cor <- matrix(runif(ncol(design_perm) * 2), nrow = ncol(design_perm))
  thinlog2 <- stats::rexp(n = n)
  design_fixed <- matrix(rep(c(0, 1), length.out = n))
  coef_fixed <- matrix(rnorm(p))
  design_obs = matrix(runif(n))

  thout <- thin_diff(mat = mat,
                     design_perm = design_perm,
                     coef_perm = coef_perm)

  summary(thout)

  thout <- thin_diff(mat = mat,
                     design_fixed = design_fixed,
                     coef_fixed = coef_fixed,
                     design_perm = design_perm,
                     coef_perm = coef_perm,
                     target_cor = target_cor,
                     design_obs = design_obs,
                     use_sva = TRUE)

  # cor(thout$sv, thout$matching_var)
  # cor(thout$sv, thout$design[, 2])
  # target_cor

  newmat <- thin_diff(mat = mat)
  expect_equal(newmat$mat, mat)
})


test_that("library thinning works", {
  set.seed(1)
  n <- 10
  p <- 1000
  mat <- matrix(1000, ncol = n, nrow = p)
  thinlog2 <- rexp(n = n, rate = 1)

  thout <- thin_lib(mat = mat, thinlog2 = thinlog2)

  empvec <- colMeans(thout$mat) / 1000
  propvec <- 2 ^ (-thinlog2)

  expect_equal(empvec, propvec, tol = 0.01)
})



test_that("fix_cor fixes the correlation", {
  set.seed(1)
  n <- 20
  p <- 2
  k <- 4
  design_perm <- matrix(rnorm(n * p), nrow = n)
  target_cor  <- matrix(runif(k * p), nrow = p)

  new_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)

  topmat <- cbind(cor(design_perm), new_cor)
  bottommat <- cbind(t(new_cor), diag(ncol(target_cor)))
  Sigma <- rbind(topmat, bottommat)

  expect_true(eigen(Sigma, only.values = TRUE)$values[k + p] > 0)
})

test_that("rmvnorm works", {
  set.seed(1)
  p <- 5
  n <- 10000
  A <- matrix(rnorm(p * p), nrow = p)
  sigma <- crossprod(A)
  mu <- matrix(1.2, nrow = 10000, ncol = 5)

  simout <- rmvnorm(mu = mu, sigma = sigma)

  expect_equal(colMeans(simout), rep(1.2, p), tol = 0.1)
  expect_equal(cov(simout), sigma, tol = 0.1)
})


test_that("permute design approximates target correlation", {
  set.seed(1)
  n <- 1000
  p <- 2
  k <- 3
  A <- matrix(rnorm((k + p)^2), nrow = p + k)
  cormat <- cov2cor(crossprod(A))
  sigma11 <- cormat[seq_len(p), seq_len(p)]
  target_cor <- cormat[seq_len(p), (p + 1):(k + p)]
  design_perm <- rmvnorm(mu = matrix(0, nrow = n, ncol = p), sigma = sigma11)
  sv <- scale(matrix(runif(k * n), nrow = n))
  target_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)
  pout <- permute_design(design_perm = design_perm, sv = sv, target_cor = target_cor, method = "marriage")
  expect_equal(cor(pout$latent_var, sv), target_cor, tol = 0.1)

  # pout <- permute_design(design_perm = design_perm, sv = sv, target_cor = target_cor, method = "optmatch")
  # cor(pout$design_perm, sv)
  # target_cor
})

test_that("permute_design orders 0's and 1's for binary designs", {
  set.seed(1)
  n <- 100
  p <- 1
  k <- 3
  design_perm <- matrix(rep(c(0, 1), length.out = n))
  sv <- scale(matrix(runif(k * n), nrow = n))
  A <- matrix(rnorm((k + p)^2), nrow = p + k)
  cormat <- cov2cor(crossprod(A))
  target_cor <- cormat[seq_len(p), (p + 1):(k + p), drop = FALSE]
  target_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)

  pout <- permute_design(design_perm = design_perm, sv = sv, target_cor = target_cor, method = "optmatch")
  expect_true(min(pout$latent_var[pout$design_perm == 1]) >= max(pout$latent_var[pout$design_perm == 0]))
})

test_that("thin_2group doesn't alter zero coef genes", {
  n <- 10
  p <- 100
  Z <- rnorm(n)
  alpha <- rnorm(p)
  mat <- round(2^(alpha %*% t(Z) + matrix(rnorm(n * p), nrow = p, ncol = n)))

  thout <- thin_2group(mat = mat)
  expect_equal(thout$mat, mat)

  thout <- thin_2group(mat = mat, prop_null = 0.5)
  expect_equal(thout$mat[abs(thout$coefmat) < 10^-6, ], mat[abs(thout$coefmat) < 10^-6, ])
})


