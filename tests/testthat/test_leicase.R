context("Lei")

test_that("Lei's case works", {
  count <- matrix(rpois(100, 10), ncol = 10)
  expect_error(poisthin(count, ngene = 10, gselect = "random"), NA)
})
