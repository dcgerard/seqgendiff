############################
## Functions relating to surrogate variables and permuting the design
## matrix in thin_diff().
############################

#' Estimate the surrogate variables.
#'
#' This will use either \code{\link[sva]{sva}} or an SVD on the residuals
#' of a regression of \code{mat} on \code{design_obs} to estimate the
#' surrogate variables.
#'
#' @inheritParams thin_diff
#' @param n_sv The number of surrogate variables.
#'
#' @author David Gerard
est_sv <- function(mat, n_sv, design_obs, use_sva = FALSE) {
  assertthat::is.count(n_sv)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.matrix(design_obs))
  assertthat::assert_that(is.logical(use_sva))
  assertthat::are_equal(ncol(mat), nrow(design_obs))
  assertthat::are_equal(length(use_sva), 1)

  matlog2 <- log2(mat + 0.5)
  Xfixed  <- cbind(design_obs, 1)
  if (use_sva & ncol(design_obs) > 0) {
    utils::capture.output(sv <- sva::sva(dat = matlog2, mod = Xfixed, n.sv = n_sv)$sv)
  } else {
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
#'   \item{Hansen, Ben B., and Stephanie Olsen Klopfer. "Optimal full matching and related designs via network flows." Journal of computational and Graphical Statistics 15, no. 3 (2006): 609-627.}
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
    design_perm <- design_perm[matchout$proposals, , drop = FALSE]
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
#' Let \eqn{W} = \code{cor(design_perm)}. Let \eqn{R} = \code{target_cor}.
#' Then the overall correlation matrix is:
#' \deqn{ \left(
#' \begin{array}{cc}
#' W  & R\\
#' R' & I_K
#' \end{array}
#' \right).
#' }
#' This function applies a multiplicative scaling factor to \eqn{R} until
#' the above matrix is positive semi-definite. That is, it finds \eqn{a}
#' between 0 and 1 such that
#' \deqn{ \left(
#' \begin{array}{cc}
#' W  & aR\\
#' aR' & I_K
#' \end{array}
#' \right)
#' }
#' is positive semi-definite.
#'
#' @inheritParams thin_diff
#' @param num_steps The number of steps between 0 and 1 to take in the
#'     grid search for the shrinkage factor. The step-size would be
#'     \code{1 / (num_steps - 1)}.
#'
#' @author David Gerard
#'
#' @examples
#' n <- 10
#' design_perm <- matrix(rep(c(0, 1), length.out = n))
#' target_cor <- matrix(seq(1, 0, length.out = 10), nrow = 1)
#' new_cor <- seqgendiff:::fix_cor(design_perm = design_perm, target_cor = target_cor)
#' new_cor / target_cor
#'
#' ## In the case of one observed covariate, the requirement is just that
#' ## the sum of squared correlations is less than or equal to one.
#' sum(target_cor ^ 2)
#' sum(new_cor ^ 2)
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

# Generate multivariate normal random samples.
#
# @param mu A matrix of means. The rows index the independent samples, the columns
#    index the variables.
# @param sigma A covariance matrix of the columns.
rmvnorm <- function(mu, sigma) {
  stopifnot(nrow(sigma) == ncol(mu))
  cholout <- chol(sigma)
  simout <- matrix(stats::rnorm(n = prod(dim(mu))), nrow = nrow(mu)) %*% cholout + mu
  return(simout)
}
