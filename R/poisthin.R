
#' Apply Poisson thinning to a matrix of count data.
#'
#' @param mat A matrix of count data. The rows index the individuals and
#'     the columns index the genes.
#' @param nsamp The number of samples to select from \code{mat}.
#' @param ngene The number of genes to select from \code{mat}.
#' @param gselect How should we select the subset of genes? Should we choose
#'     the \code{ngene} most expressed genes (\code{"max"}), a random sample
#'     of the genes (\code{"random"}), a random sample of the most expressed
#'     genes (\code{"rand_max"}), or a user-provided list (\code{"custom"})?
#'     If \code{"custom"}, then \code{gvec} should be specified. Expression levels
#'     of a gene are measured by median expression across individuals with ties broken
#'     by mean expression.
#' @param gvec A logical of length \code{ncol(mat)}. A \code{TRUE} in position \eqn{i}
#'     indicates inclusion into the smaller dataset. Hence, \code{sum(gvec)} should
#'     equal \code{ngene}.
#' @param signal_fun A function that returns the signal. This should take as
#'     input \code{n} for the number of samples to return and then return only
#'     a vector of samples.
#' @param signal_params A list of additional arguments to pass to \code{signal_fun}.
#' @param skip_gene The number of maximally expressed genes to skip.
#' @param prop_null The proportion of genes that are null.
#'
#' @return A list with the following elements:
#'
#' @author David Gerard
#'
#' @export
poisthin <- function(mat, nsamp = nrow(mat), ngene = ncol(mat),
                     gselect = c("max", "random", "rand_max", "custom"),
                     gvec = NULL,
                     skip_gene = 0,
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 1),
                     prop_null = 1) {

  ## Check Input -------------------------------------------------------------
  nsamp     <- as.integer(nsamp)
  ngene     <- as.integer(ngene)
  skip_gene <- as.integer(skip_gene)

  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nsamp <= nrow(mat))
  assertthat::assert_that(ngene + skip_gene <= ncol(mat))
  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::assert_that(prop_null >= 0, prop_null <= 1)

  gselect <- match.arg(gselect)

  if (gselect == "custom") {
    assertthat::assert_that(is.logical(gvec))
    assertthat::are_equal(length(gvec), ncol(mat))
    assertthat::are_equal(sum(gvec), ngene)
  } else {
    if (!is.null(gvec)) {
      warning('gvec is specified but being ignored since gselect is not "custom"')
    }
  }

  ## subset matrix -----------------------------------------------------------
  if (gselect == "max" | gselect == "rand_max") {
    med_express <- apply(mat, 2, stats::median)
    mean_express <- colMeans(mat)
    order_vec <- order(med_express, mean_express)
  }

  if (gselect == "max") {
    gindices <- order_vec[1:ngene]
  } else if (gselect == "rand_max") {
    first_zero <- match(0, med_express)
    max_gene <- min(c(first_zero, 2 * ngene, ncol(mat)), na.rm = TRUE)
    gindices <- sample(x = order_vec[1:max_gene], size = ngene)
  } else if (gselect == "random") {
    gindices <- sample(x = sample(1:ncol(mat)), size = ngene)
  } else if (gselect == "custom") {
    gindices <- (1:ncol(mat))[gvec]
  }

  samp_indices <- sample(1:nrow(mat), size = nsamp)

  submat <- mat[samp_indices, gindices, drop = FALSE]
  group_indicator <- rep(0, length = nsamp)
  group_indicator[sample(1:nsamp, size = floor(nsamp / 2))] <- 1

  ## Draw signal -------------------------------------------------------------
  nsignal <- round(ngene * (1 - prop_null))
  if (nsignal > 0) {
    signal_params$n <- nsignal
    signal_vec      <- do.call(what = signal_fun, args = signal_params) ## log2-fold change
    which_signal    <- sample(1:ncol(submat), nsignal) # location of signal
    sign_vec        <- sign(signal_vec) # sign of signal
    bin_probs       <- 2 ^ -abs(signal_vec) # binomial prob
    for (gn in 1:length(signal_vec)) {
      if (sign_vec[gn] == 1) {
        current_count <- submat[group_indicator == 1, which_signal[gn]]
        submat[group_indicator == 1, which_signal[gn]] <-
          sapply(current_count, FUN = stats::rbinom, n = 1, prob = bin_probs[gn])
      } else if (sign_vec[gn] == -1) {
        current_count <- submat[group_indicator == 0, which_signal[gn]]
        submat[group_indicator == 0, which_signal[gn]] <-
          sapply(current_count, FUN = stats::rbinom, n = 1, prob = bin_probs[gn])
      }
    }
    beta <- rep(0, ngene)
    beta[which_signal] <- signal_vec
  } else if (nsignal == 0 & abs(prop_null - 1) < 10 ^ -6) {
    warning('no genes were given signal since (1 - prop_null) * ngene was very close to zero')
  } else {
    beta <- rep(0, ngene)
  }

  X <- stats::model.matrix(~group_indicator)
  return_list <- list(Y = submat, X = X, beta = beta)

  return(return_list)
}
