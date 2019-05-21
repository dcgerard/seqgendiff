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

  thout <- thin_diff(mat = mat,
                     design_perm = design_perm,
                     coef_perm = coef_perm)

  thout <- thin_diff(mat = mat,
                     design_fixed = design_fixed,
                     coef_fixed = coef_fixed,
                     design_perm = design_perm,
                     coef_perm = coef_perm,
                     target_cor = target_cor,
                     use_sva = TRUE)

  # cor(thout$sv, thout$matching_var)
  # cor(thout$sv, thout$design[, 2])
  # target_cor

  newmat <- thin_diff(mat = mat)
  expect_equal(newmat$mat, mat)
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


  # n <- 10
  # design_perm <- rmvnorm(mu = matrix(0, nrow = n, ncol = p), sigma = sigma11)
  # sv <- scale(matrix(runif(k * n), nrow = n))
  # target_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)
  #
  # itermax <- 1000
  # corarray <- array(0, dim = c(p, k, itermax))
  # for (index in seq_len(itermax)) {
  #   pout <- permute_design(design_perm = design_perm, sv = sv, target_cor = target_cor, method = "optmatch")
  #   corarray[,,index] <- cor(pout$design_perm, sv)
  # }
  # apply(corarray, c(1, 2), median)
  # target_cor
})






