#########################
## General Thinning functions
#########################

#' Base Poisson thinning function.
#'
#' Given a matrix of counts (\eqn{Y}), a design matrix (\eqn{X}),
#' and a matrix of coefficients (\eqn{B}), \code{thin_diff} will generate a new
#' matrix of counts such that \eqn{E[log_2(Y)] \approx BX' + u1'}, where
#' \eqn{u} is some vector of intercept coefficients. This
#' function is used by all other thinning functions.
#'
#' @inheritParams thin_diff
#' @param designmat A design matrix. The rows index the samples and the columns
#'    index the variables. The intercept should \emph{not} be included.
#' @param coefmat A matrix of coefficients. The rows index the genes and the
#'    columns index the samples.
#'
#' @export
#'
#' @author David Gerard
#'
#' @examples
#' ## Simulate data from given matrix of counts
#' set.seed(1)
#' nsamp <- 10
#' ngene <- 1000
#' Y <- matrix(stats::rpois(nsamp * ngene, lambda = 100), nrow = ngene)
#' X <- matrix(rep(c(0, 1), length.out = nsamp))
#' B <- matrix(seq(3, 0, length.out = ngene))
#' Ynew <- thin_base(mat = Y, designmat = X, coefmat = B)
#'
#' ## Demonstrate how the log2 effect size is B
#' Bhat <- coefficients(lm(t(log2(Ynew)) ~ X))["X", ]
#' plot(Bhat, B)
#' abline(0, 1, col = 2, lwd = 2)
thin_base <- function(mat, designmat, coefmat, relative = TRUE) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.matrix(designmat))
  assertthat::assert_that(is.matrix(coefmat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::assert_that(is.numeric(designmat))
  assertthat::assert_that(is.numeric(coefmat))
  assertthat::are_equal(nrow(mat), nrow(coefmat))
  assertthat::are_equal(ncol(mat), nrow(designmat))
  assertthat::are_equal(ncol(designmat), ncol(coefmat))
  assertthat::assert_that(is.logical(relative))
  assertthat::are_equal(1L, length(relative))
  stopifnot(mat >= 0)

  ## Thin ---------------------------------------------------------------------
  meanmat <- tcrossprod(coefmat, designmat)
  maxvec  <- apply(meanmat, 1, max)
  if (!relative) {
    if (any(maxvec > 0)) {
      stop("thin_base: tcrossprod(coefmat, designmat) produced positive entries.")
    }
    qmat <- 2 ^ meanmat
  } else {
    qmat <- 2 ^ (meanmat - maxvec)
  }
  newmat      <- stats::rbinom(n = prod(dim(mat)), size = mat, prob = qmat)
  dim(newmat) <- dim(mat)
  return(newmat)
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
thin_lib <- function(mat, thinlog2, relative = FALSE) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::are_equal(length(thinlog2), ncol(mat))
  thinlog2 <- c(thinlog2)
  assertthat::assert_that(is.numeric(thinlog2))
  stopifnot(thinlog2 >= 0)

  thout <- thin_diff(mat          = mat,
                     design_fixed = matrix(-thinlog2, ncol = 1),
                     coef_fixed   = matrix(1, nrow = nrow(mat), ncol = 1),
                     relative     = relative)

  return(thout)
}

#' Poisson thinning in the two-group model.
#'
#'
#'
#' @inheritParams thin_diff
#' @param prop_null The proportion of genes that are null.
#' @param signal_fun A function that returns the signal. This should take as
#'     input \code{n} for the number of samples to return and then return only
#'     a vector of samples. Additional parameters may be passed through
#'     \code{signal_params}.
#' @param signal_params A list of additional arguments to pass to
#'     \code{signal_fun}.
#' @param group_prop The proportion of individuals that are in group 1.
#' @param corvec A vector of target correlations. \code{corvec[i]} is the
#'     target correlation of the latent group assignment vector with the
#'     ith latent confounder. The default is to set this to \code{NULL},
#'     in which case group assignment is made independently of any
#'     unobserved confounding.
#'
#' @author David Gerard
thin_2group <- function(mat,
                        prop_null     = 1,
                        signal_fun    = stats::rnorm,
                        signal_params = list(mean = 0, sd = 1),
                        group_prop    = 0.5,
                        corvec        = NULL) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::are_equal(1L, length(prop_null), length(signal_fun), length(group_prop))
  assertthat::assert_that(prop_null <= 1, prop_null >= 0, group_prop <= 1, group_prop >= 0)
  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::assert_that(is.null(signal_params$n))
  ngene <- nrow(mat)
  nsamp <- ncol(mat)


  ## Generate coef ------------------------------------------------------------
  coef_perm <- matrix(0, nrow = ngene, ncol = 1)
  signal_params$n <- round(ngene * (1 - prop_null))
  if (signal_params$n > 0) {
    signal_vec <- c(do.call(what = signal_fun, args = signal_params))
    coef_perm[sample(seq_len(ngene), size = signal_params$n), 1] <- signal_vec
  }

  ## Generate design ----------------------------------------------------------
  numtreat <- round(group_prop * nsamp)
  if (numtreat == 0) {
    design_perm <- matrix(0, ncol = 1, nrow = nsamp)
  } else if (numtreat == nsamp) {
    design_perm <- matrix(1, ncol = 1, nrow = nsamp)
  } else {
    design_perm <- matrix(0, ncol = 1, nrow = nsamp)
    design_perm[sample(seq_len(nsamp), size = numtreat), ] <- 1
  }

  ## Generate target correlation ----------------------------------------------
  if (!is.null(corvec)) {
    target_cor <- matrix(corvec, nrow = 1)
  } else {
    target_cor <- NULL
  }

  ## Thin ---------------------------------------------------------------------
  thout <- thin_diff(mat         = mat,
                     design_perm = design_perm,
                     coef_perm   = coef_perm,
                     target_cor  = target_cor)

  return(thout)
}

#' Poisson thinning for differential expression analysis.
#'
#' @param mat A matrix of counts. The rows index the genes and the columns
#'     index the samples (as is usual in RNA-seq).
#' @param design_fixed A design matrix whose rows are fixed and not permuted.
#'     The rows index the samples and the columns index the variables.
#'     The intercept should \emph{not} be included.
#' @param coef_fixed The coefficients corresponding to \code{design_fixed}.
#'     The rows index the genes and the columns index the variables.
#' @param design_perm A design matrix whose rows are to be permuted (thus
#'     controlling the amount by which they are correlated with the confounders).
#'     The rows index the samples and the columns index the variables.
#'     The intercept should \emph{not} be included.
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
#'     The intercept should \emph{not} be included.
#' @param relative A logical. Should we apply relative thinning (\code{TRUE})
#'     or absolute thinning (\code{FALSE}). Only experts should change
#'     the default.
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
                      design_obs   = NULL,
                      relative     = TRUE) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::assert_that(is.logical(use_sva))
  assertthat::assert_that(is.logical(relative))
  assertthat::are_equal(1L, length(use_sva), length(relative))
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
    new_cor <- fix_cor(design_perm = design_perm,
                       target_cor  = target_cor)

    ## Estimate hidden confounders ----------------------
    n_sv <- ncol(new_cor)
    sv <- est_sv(mat          = mat,
                 n_sv         = n_sv,
                 design_fixed = design_fixed,
                 design_obs   = design_obs,
                 use_sva      = use_sva)

    ## Permute design matrix ----------------------------
    pout <- permute_design(design_perm = design_perm,
                           sv          = sv,
                           target_cor  = new_cor)
    design_perm <- pout$design_perm
    latent_var  <- pout$latent_var
  }

  ## Make overall design and coef ---------------------------------------------
  designmat        <- cbind(design_fixed, design_perm)
  class(designmat) <- "numeric"
  coefmat          <- cbind(coef_fixed, coef_perm)
  class(coefmat)   <- "numeric"

  ## Thin ---------------------------------------------------------------------
  newmat <- thin_base(mat       = mat,
                      designmat = designmat,
                      coefmat   = coefmat,
                      relative  = relative)

  retval <- list(mat          = newmat,
                 designmat    = designmat,
                 coefmat      = coefmat,
                 sv           = sv,
                 cor          = new_cor,
                 matching_var = latent_var)

  class(retval) <- "ThinData"

  return(retval)
}

#' Estimate the surrogate variables.
#'
#' @inheritParams thin_diff
#' @param n_sv The number of surrogate variables.
#'
#' @author David Gerard
est_sv <- function(mat, n_sv, design_fixed, design_obs, use_sva = FALSE) {
  assertthat::is.count(n_sv)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.matrix(design_fixed))
  assertthat::assert_that(is.matrix(design_obs))
  assertthat::assert_that(is.logical(use_sva))
  assertthat::are_equal(ncol(mat), nrow(design_fixed), nrow(design_obs))
  assertthat::are_equal(length(use_sva), 1)

  matlog2 <- log2(mat + 0.5)
  if (use_sva & ncol(design_fixed) > 0) {
    utils::capture.output(sv <- sva::sva(dat = matlog2, mod = cbind(design_fixed, design_obs, 1), n.sv = n_sv)$sv)
  } else {
    Xfixed <- cbind(design_fixed, design_obs, 1)
    sv <- svd(matlog2 %*% (diag(nrow(Xfixed)) - Xfixed %*% solve(t(Xfixed) %*% Xfixed) %*% t(Xfixed)), nv = n_sv, nu = 0)$v
  }
  sv <- sv * sqrt(nrow(sv - 1))
  return(sv)
}

#' Permute the design matrix so that it is approximately correlated with
#' the surrogate variables.
#'
#' @inheritParams thin_diff
#' @param sv A matrix of surrogate variables
#' @param method Should we use the optimal matching technique from Hansen and
#'     Klopfer (2006) (\code{"optmatch"}) or the Gale-Shapley algorithm
#'     for stable marriages (\code{marriage}) (Gale and Shapley, 1962)
#'     as implemented in the matchingR package.
#'     The \code{"optmatch"} method works almost uniformly better in practice,
#'     but does take a lot more computational time if you have, say, 1000
#'     samples.
#'
#' @references
#' \itemize{
#'   \item{Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and related designs via network flows, JCGS 15 609-627.}
#'   \item{Gale, David, and Lloyd S. Shapley. "College admissions and the stability of marriage." The American Mathematical Monthly 69, no. 1 (1962): 9-15.}
#' }
#'
#' @author David Gerard
permute_design <- function(design_perm, sv, target_cor, method = c("optmatch", "marriage")) {
  ## Check input --------------------------------------------------------------
  assertthat::are_equal(nrow(design_perm), nrow(sv))
  assertthat::are_equal(ncol(design_perm), nrow(target_cor))
  assertthat::are_equal(ncol(sv), ncol(target_cor))
  method <- match.arg(method)
  nsamp <- nrow(design_perm)
  ## in case a non-standard sv is given ---------------------------------------
  sv <- scale(sv)

  ## Generate latent factors --------------------------------------------------
  sigma11 <- stats::cor(design_perm)
  sigma12 <- target_cor
  sigma_cond <- sigma11 - sigma12 %*% t(sigma12)
  mu_cond <- sv %*% t(sigma12)
  latent_var <- rmvnorm(mu = mu_cond, sigma = sigma_cond)

  ## Get permutations ---------------------------------------------------------
  distmat <- as.matrix(pdist::pdist(X = scale(design_perm), Y = scale(latent_var)))
  if (method == "optmatch") {
    dimnames(distmat) <- list(treated = paste0("O", seq_len(nsamp)), control = paste0("L", seq_len(nsamp)))
    suppressWarnings(matchout <- optmatch::pairmatch(distmat))
    ogroup <- matchout[attributes(matchout)$contrast.group]
    lgroup <- matchout[!attributes(matchout)$contrast.group]
    design_perm <- design_perm[match(lgroup, ogroup), , drop = FALSE]
  } else if (method == "marriage") {
    matchout <- matchingR::galeShapley.marriageMarket(proposerUtils = -1 * t(distmat), reviewerUtils = -1 * distmat)
    design_perm <- design_perm[matchout$proposals, ]
  }
  return(list(design_perm = design_perm, latent_var = latent_var))
}

#' Estimate the effective correlation
#'
#' Will return the actual correlation between the design matrix and the
#' surrogate variables when you use \code{\link{permute_design}}.
#'
#' @inheritParams thin_diff
#' @inheritParams permute_design
#' @param iternum The total number of simulated correlations to consider.
#'
#' @author David Gerard
effective_cor <- function(design_perm, sv, target_cor, method = c("optmatch", "marriage"), iternum = 1000) {
  method <- match.arg(method)
  target_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)
  itermax <- 1000
  corarray <- array(0, dim = c(ncol(design_perm), ncol(sv), itermax))
  for (index in seq_len(itermax)) {
    pout <- permute_design(design_perm = design_perm, sv = sv, target_cor = target_cor, method = method)
    corarray[,,index] <- stats::cor(pout$design_perm, sv)
  }
  truecor <- apply(corarray, c(1, 2), mean)
  return(truecor)
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

