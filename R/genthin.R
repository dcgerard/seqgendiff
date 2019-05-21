#########################
## General Thinning functions
#########################



#' Basic Poisson thinning function.
#'
#' Given a matrix of counts, a design matrix, and a matrix
#'
#' @inheritParams thin_diff
#' @param design A design matrix. The rows index the samples and the columns
#'    index the variables.
#' @param coef A matrix of coefficients. The rows index the genes and the
#'    columns index the samples.
#'
#' @author David Gerard
thin_base <- function(mat, design, coef) {
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.matrix(design))
  assertthat::assert_that(is.matrix(coef))
  assertthat::assert_that(is.numeric(mat))
  assertthat::assert_that(is.numeric(design))
  assertthat::assert_that(is.numeric(coef))
  assertthat::are_equal(nrow(mat), nrow(coef))
  assertthat::are_equal(ncol(mat), nrow(design))
  assertthat::are_equal(ncol(design), ncol(coef))
  stopifnot(mat >= 0)

  meanmat <- tcrossprod(coef, design)
  qmat <- 2 ^ (meanmat - apply(meanmat, 1, max))
  newmat <- stats::rbinom(n = prod(dim(mat)), size = mat, prob = qmat)
  dim(newmat) <- dim(mat)
  return(newmat)
}

#' Poisson thinning for differential expression analysis.
#'
#' @param mat A matrix of counts. The rows index the genes and the columns
#'     index the samples (as is usual in RNA-seq).
#' @param design_fixed A design matrix whose rows are fixed and not permuted.
#'     The rows index the samples and the columns index the variables.
#' @param coef_fixed The coefficients corresponding to \code{design_fixed}.
#'     The rows index the genes and the columns index the variables.
#' @param design_perm A design matrix whose rows are to be permuted (thus
#'     controlling the amount by which they are correlated with the confounders).
#'     The rows index the samples and the columns index the variables.
#' @param coef_perm The coefficients corresponding to \code{design_perm}.
#'     The rows index the genes and the columns index the variables.
#' @param target_cor A matrix of target correlation betweens the variables in
#'     \code{design_perm} and the unobserved confounders. The rows index the
#'     observed covariates and the columns index the unobserved confounders.
#'     The number of columns indicates the number of hidden confounders.
#'     Set this to \code{NULL} to indicate uncorrelated.
#' @param use_sva A logical. Should we use SVA using \code{design_fixed}
#'     and \code{design_obs} to estimate the hidden covariates (\code{TRUE})
#'     or should we just do an SVD on \code{log2(mat + 0.5)} after
#'     subtracing row means (\code{FALSE})? Defaults to \code{FALSE}.
#' @param design_obs A matrix of observed covariates that we are NOT to be a
#'     part of the signal generating process. Only used in estimating the
#'     confounders if \code{use_sva = TRUE}.
#'
#' @export
#'
#' @author David Gerard
thin_diff <- function(mat,
                      design_fixed = NULL,
                      coef_fixed   = NULL,
                      design_perm  = NULL,
                      coef_perm    = NULL,
                      target_cor   = NULL,
                      use_sva      = FALSE,
                      design_obs   = NULL) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  stopifnot(mat >= 0)
  assertthat::assert_that(is.logical(use_sva))
  ngene <- nrow(mat)
  nsamp <- ncol(mat)

  if (is.null(design_fixed) | is.null(coef_fixed)) {
    design_fixed <- matrix(ncol = 0L, nrow = nsamp)
    class(design_fixed) <- "numeric"
    coef_fixed   <- matrix(ncol = 0L, nrow = ngene)
    class(coef_fixed) <- "numeric"
  }
  if (is.null(design_perm) | is.null(coef_perm)) {
    design_perm <- matrix(ncol = 0L, nrow = nsamp)
    class(design_perm) <- "numeric"
    coef_perm   <- matrix(ncol = 0L, nrow = ngene)
    class(coef_perm) <- "numeric"
  }
  if (is.null(design_obs)) {
    design_obs <- matrix(ncol = 0L, nrow = nsamp)
    class(design_obs) <- "numeric"
  }

  assertthat::are_equal(ncol(mat), nrow(design_fixed), nrow(design_perm), nrow(design_obs))
  assertthat::are_equal(nrow(mat), nrow(coef_fixed), nrow(coef_perm))
  assertthat::are_equal(ncol(design_fixed), ncol(coef_fixed))
  assertthat::are_equal(ncol(design_perm), ncol(coef_perm))

  ## Permute ------------------------------------------------------------------
  if (is.null(target_cor) | nrow(design_perm) == 0) {
    design_perm <- design_perm[sample(seq_len(nrow(design_perm))), ]
    new_cor <- NULL
    latent_var <- matrix(ncol = 0L, nrow = nsamp)
    class(latent_var) <- "numeric"
    sv <- matrix(ncol = 0L, nrow = nsamp)
    class(sv) <- "numeric"
  } else {
    ## Fix target correlation ---------------------------
    new_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)

    ## Estimate hidden confounders ----------------------
    n_sv <- ncol(new_cor)
    matlog2 <- log2(mat + 0.5)
    if (use_sva & ncol(design_fixed) > 0) {
      utils::capture.output(sv <- sva::sva(dat = matlog2, mod = cbind(design_fixed, design_obs, 1), n.sv = n_sv)$sv)
    } else {
      Xfixed <- cbind(design_fixed, 1)
      sv <- svd(matlog2 %*% (diag(nrow(Xfixed)) - Xfixed %*% solve(t(Xfixed) %*% Xfixed) %*% t(Xfixed)), nv = 2, nu = 0)$v
    }
    sv <- sv * sqrt(nrow(sv))

    ## Generate latent factors --------------------------
    sigma11 <- stats::cor(design_perm)
    sigma12 <- target_cor
    sigma_cond <- sigma11 - sigma12 %*% t(sigma12)
    mu_cond <- sv %*% t(sigma12)
    latent_var <- rmvnorm(mu = mu_cond, sigma = sigma_cond)

    ## Get permutations ---------------------------------
    distmat <- as.matrix(pdist::pdist(X = design_perm, Y = latent_var))
    dimnames(distmat) <- list(treated = paste0("O", seq_len(nsamp)), control = paste0("L", seq_len(nsamp)))
    suppressWarnings(matchout <- optmatch::pairmatch(distmat))
    ogroup <- matchout[attributes(matchout)$contrast.group]
    lgroup <- matchout[!attributes(matchout)$contrast.group]
    design_perm <- design_perm[match(lgroup, ogroup), , drop = FALSE]
  }

  ## Make overall design and coef ---------------------------------------------
  design <- cbind(design_fixed, design_perm)
  class(design) <- "numeric"
  coef   <- cbind(coef_fixed, coef_perm)
  class(coef) <- "numeric"

  ## Thin ---------------------------------------------------------------------
  newmat <- thin_base(mat    = mat,
                      design = design,
                      coef   = coef)

  return(list(mat = newmat, design = design, coef = coef, sv = sv, cor = new_cor, matching_var = latent_var))
}

#' Poisson thinning for altering library size.
#'
#' @inheritParams thin_diff
#' @param thinlog2 A vector. Element i is the amount to thin (on the log2 scale). For
#'     example, a value of 0 means that we do not thin, a value of 1 means
#'     that we thin by a factor of 2, a value of 2 means we thin by a factor
#'     of 4, etc.
#'
#' @export
#'
#' @author David Gerard
thin_lib <- function(mat, thinlog2) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  stopifnot(mat >= 0)
  assertthat::are_equal(length(thinlog2), ncol(mat))
  thinlog2 <- c(thinlog2)
  assertthat::assert_that(is.numeric(thinlog2))
  stopifnot(thinlog2 >= 0)

  newmat <- thin_base(mat    = mat,
                      design = matrix(-thinlog2, ncol = 1),
                      coef   = matrix(1, nrow = nrow(mat), ncol = 1))

  return(newmat)
}



#' Fixes an invalid target correlation.
#'
#' Shrinks the target correlation using a uniform scaling facter so that
#' the overall correlation matrix is positive semi-definite.
#'
#' @inheritParams thin_diff
#' @param num_steps The number of steps between 0 and 1 to take in the
#'     grid search for the shrinkage factor. The step-size would be
#'     \code{1 / (num_steps - 1)}.
#'
#' @author David Gerard
fix_cor <- function(design_perm, target_cor, num_steps = 51) {
  ## Check input --------------------------------------------------------------
  stopifnot(is.matrix(design_perm))
  stopifnot(is.matrix(target_cor))
  stopifnot(ncol(design_perm) == nrow(target_cor))
  stopifnot(target_cor >= -1, target_cor <= 1)
  assertthat::is.count(num_steps)
  stopifnot(num_steps > 1)

  p <- ncol(design_perm)
  ## Calculate top left correlation -------------------------------------------
  top_left_cor <- stats::cor(design_perm)

  ## Calculate RR^T -----------------------------------------------------------
  RRt <- tcrossprod(target_cor)

  ## Test if already psd ------------------------------------------------------
  eout <- eigen(top_left_cor - RRt, only.values = TRUE)
  if (eout$values[p] >= 0) {
    return(target_cor)
  }

  shrink_vec <- seq(1, 0, length.out = num_steps)
  valid <- FALSE
  step_index <- 1
  while(!valid) {
    step_index <- step_index + 1
    current_shrink <- shrink_vec[step_index]
    eout <- eigen(top_left_cor - current_shrink * RRt, only.values = TRUE)
    if (eout$values[p] >= 0) {
      valid <- TRUE
    }
  }

  return(sqrt(current_shrink) * target_cor)
}

# Generate multivariate normal random variable.
rmvnorm <- function(mu, sigma) {
  stopifnot(nrow(sigma) == ncol(mu))
  cholout <- chol(sigma)
  simout <- matrix(stats::rnorm(n = prod(dim(mu))), nrow = nrow(mu)) %*% cholout + mu
  return(simout)
}

