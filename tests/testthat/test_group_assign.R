context("Group Assignment")

test_that("length of corvec and nfac are 0", {
  set.seed(2)
  n <- 10
  p <- 100
  Z <- rnorm(n)
  alpha <- rnorm(p)
  Y <- round(2 ^ (Z %*% t(alpha) + matrix(rnorm(n * p), nrow = n, ncol = p)))

  expect_error(corassign(mat = Y, nfac = 0, corvec = 0.2))

  expect_message(corassign(mat = Y, corvec = c(0.2, 0.2, 0.2)))

  expect_message(cout <- corassign(mat = Y, corvec = 0, return = "full"))
  expect_equal(cout$nfac, ncol(cout$facmat), length(cout$corvec))

  cout <- corassign(mat = Y, nfac = 0, return = "full")
  expect_equal(cout$nfac, ncol(cout$facmat), length(cout$corvec))

  expect_error(corassign(mat = Y, nfac = 1, corvec = 1.1))

  cout <- corassign(mat = Y, nfac = 2, corvec = c(0.1, 0.2), return = "full")
  expect_equal(cout$nfac, ncol(cout$facmat), length(cout$corvec))

  gout <- corassign(mat = Y, nfac = 2, corvec = c(0, 0))
  gout <- corassign(mat = Y)
  gout <- corassign(mat = Y, nfac = 2)
  expect_message(gout <- corassign(mat = Y, corvec = 0.1))
})

test_that("uncorassign works", {
  set.seed(1)
  fout1 <- uncorassign(n = 100000, return = "group")
  fout2 <- uncorassign(n = 100000, return = "full")

  expect_equal(mean(fout1$x), 0.5, tolerance = 0.01)
  expect_equal(mean(fout2$x), 0.5, tolerance = 0.01)
})


test_that("getting correlation close to what you expect", {
  set.seed(1)
  n <- 10000
  p <- 10
  Z <- rnorm(n)
  alpha <- rnorm(p)
  Y <- round(2 ^ (Z %*% t(alpha) + matrix(rnorm(n * p), nrow = n, ncol = p)))
  tol <- 0.02

  cout5 <- corassign(mat = Y, nfac = 1, corvec = 0.5, return = "full")
  expect_equal(cor(cout5$facmat[, 1, drop = TRUE], cout5$groupfac), 0.5,
               tolerance = tol)

  cout25 <- corassign(mat = Y, nfac = 1, corvec = 0.25, return = "full")
  expect_equal(cor(cout25$facmat[, 1, drop = TRUE], cout25$groupfac), 0.25,
               tolerance = tol)

  cout0 <- corassign(mat = Y, nfac = 1, corvec = 0, return = "full")
  expect_equal(cor(cout0$facmat[, 1, drop = TRUE], cout0$groupfac), 0,
               tolerance = tol)

  coutn5 <- corassign(mat = Y, nfac = 2, corvec = c(-0.5, 0.1), return = "full")
  expect_equal(cor(coutn5$facmat, coutn5$groupfac),
               matrix(c(-0.5, 0.1), ncol = 1),
               tolerance = tol)
})



test_that("get about the right proportions in 'frac' and 'random'", {
  set.seed(1)
  n <- 10000
  p <- 10
  Z <- rnorm(n)
  alpha <- rnorm(p)
  Y <- round(2 ^ (Z %*% t(alpha) + matrix(rnorm(n * p), nrow = n, ncol = p)))
  pout <- poisthin(mat = Y, group_assign = "frac", group_prop = 0.2)
  expect_equal(sum(pout$X[, 2]), 2000)

  pout <- poisthin(mat = Y, group_assign = "random", group_prop = 0.2)
  expect_equal(mean(pout$X[, 2]), 0.2, tol = 0.01)
})
