context("select_gene")

test_that("select_counts works", {
  set.seed(2)
  n <- 100
  p <- 1000
  Z <- rnorm(n)
  alpha <- rnorm(p)
  mat <- rpois(n = n * p, lambda = t(2 ^ (Z %*% t(alpha) + matrix(rnorm(n * p), nrow = n, ncol = p)) * 100))
  dim(mat) <- c(p, n)

  nsamp        <- 10
  ngene        <- 100
  gselect      <- "max"
  gvec         <- c(rep(TRUE, ngene), rep(FALSE, nrow(mat) - ngene))
  filter_first <- TRUE
  nskip        <- 1

  meanvec <- rowMeans(mat)
  medvec  <- apply(mat, 1, median)

  ## Test mean_max ------------------------------------------------------------
  submat <- select_counts(mat = mat,
                          nsamp = nsamp,
                          ngene = ngene,
                          gselect = "mean_max",
                          filter_first = FALSE,
                          nskip = 20)
  which_gene <- as.numeric(rownames(submat))
  ordermean <- order(meanvec, decreasing = TRUE)
  ordermed  <- order(medvec, meanvec, decreasing = TRUE)
  expect_false(any(ordermed[seq_len(20)] %in% which_gene))
  expect_equal(sort(ordermean[!(ordermean %in% ordermed[seq_len(20)])][seq_len(ngene)]),
               which_gene)

  ## Test gvec ----------------------------------------------------------------
  submat <- select_counts(mat = mat,
                          nsamp = nsamp,
                          ngene = ngene,
                          gselect = "custom",
                          gvec = gvec,
                          filter_first = FALSE,
                          nskip = 0L)

  ## Test 0 selections --------------------------------------------------------
  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = ngene, gselect = "max")
  expect_equal(ncol(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = nsamp, ngene = 0L, gselect = "max")
  expect_equal(nrow(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = 0L, gselect = "max")
  expect_equal(dim(countmat), c(0L, 0L))

  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = ngene, gselect = "mean_max")
  expect_equal(ncol(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = nsamp, ngene = 0L, gselect = "mean_max")
  expect_equal(nrow(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = 0L, gselect = "mean_max")
  expect_equal(dim(countmat), c(0L, 0L))

  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = ngene, gselect = "random")
  expect_equal(ncol(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = nsamp, ngene = 0L, gselect = "random")
  expect_equal(nrow(countmat), 0L)
  countmat <- select_counts(mat = mat, nsamp = 0L, ngene = 0L, gselect = "random")
  expect_equal(dim(countmat), c(0L, 0L))

  countmat <- select_counts(mat = mat)
  expect_true(all(mat == countmat))

  ## Test errors and warnings -------------------------------------------------
  expect_warning(select_counts(mat = mat, gvec = gvec, gselect = "max", ngene = ngene))
  expect_warning(select_counts(mat = mat, gvec = gvec, gselect = "custom", nskip = 2, ngene = ngene))
  expect_error(select_counts(mat = mat, gvec = gvec, gselect = "custom", ngene = ngene, filter_first = TRUE))

  if (requireNamespace("edgeR", quietly = TRUE)) {
    countmat <- select_counts(mat = mat, filter_first = TRUE)
  }

})
