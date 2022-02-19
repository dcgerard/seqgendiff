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
#' @return A matrix of estimated surrogate variables. The columns index the
#'     surrogate variables and the rows index the individuals. The surrogate
#'     variables are centered and scaled to have mean 0 and variance 1.
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
#' @param method Should we use the Gale-Shapley algorithm
#'     for stable marriages (\code{"marriage"}) (Gale and Shapley, 1962)
#'     as implemented in the matchingR package, or the Hungarian algorithm
#'     (Papadimitriou and Steiglitz, 1982) (\code{"hungarian"})
#'     as implemented in the clue package (Hornik, 2005)? The
#'     Hungarian method almost always works better, so is the default.
#'
#' @references
#' \itemize{
#'   \item{Gale, David, and Lloyd S. Shapley. "College admissions and the stability of marriage." \emph{The American Mathematical Monthly} 69, no. 1 (1962): 9-15. \doi{10.1080/00029890.1962.11989827}.}
#'   \item{C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization: Algorithms and Complexity. Englewood Cliffs: Prentice Hall.}
#'   \item{Hornik K (2005). "A CLUE for CLUster Ensembles." \emph{Journal of Statistical Software}, 14(12). \doi{10.18637/jss.v014.i12}. \doi{10.18637/jss.v014.i12}.}
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{design_perm}}{A row-permuted version of the user-provided
#'       \code{design_perm}.}
#'   \item{\code{latent_var}}{A matrix of the latent variables on which
#'       \code{design_perm} was matched.}
#' }
#'
#' @author David Gerard
permute_design <- function(design_perm, sv, target_cor, method = c("hungarian", "marriage")) {
  ## Check input --------------------------------------------------------------
  assertthat::are_equal(nrow(design_perm), nrow(sv))
  assertthat::are_equal(ncol(design_perm), nrow(target_cor))
  assertthat::are_equal(ncol(sv), ncol(target_cor))
  nsamp <- nrow(design_perm)
  ## in case a non-standard sv is given ---------------------------------------
  sv <- scale(sv)

  method <- match.arg(method)

  ## Generate latent factors --------------------------------------------------
  sigma11 <- stats::cor(design_perm)
  sigma12 <- target_cor
  sigma_cond <- sigma11 - sigma12 %*% t(sigma12)
  mu_cond <- sv %*% t(sigma12)
  latent_var <- rmvnorm(mu = mu_cond, sigma = sigma_cond)

  ## Get permutations ---------------------------------------------------------
  distmat <- as.matrix(pdist::pdist(X = scale(design_perm), Y = scale(latent_var)))
  if (method == "marriage") {
    matchout <- matchingR::galeShapley.marriageMarket(proposerUtils = -1 * distmat, reviewerUtils = -1 * t(distmat))
    design_perm <- design_perm[matchout$proposals, , drop = FALSE]
  } else if (method == "hungarian") {
    clue_out <- clue::solve_LSAP(x = t(distmat), maximum = FALSE)
    design_perm <- design_perm[as.numeric(clue_out), , drop = FALSE]
  }
  return(list(design_perm = design_perm, latent_var = latent_var))
}

#' Estimates the effective correlation.
#'
#' Will return the estimated correlation between the design matrix and the
#' surrogate variables when you assign a target correlation. The method is
#' described in detail in Gerard (2020).
#'
#' This function permutes the rows of \code{design_perm} many times, each
#' time calculating the Pearson correlation between the columns of
#' \code{design_perm} and the columns of \code{sv}. It then returns the
#' averages of these Pearson correlations. The permutation is done
#' using \code{\link{permute_design}}.
#'
#' @inheritParams thin_diff
#' @inheritParams permute_design
#' @param iternum The total number of simulated correlations to consider.
#' @param calc_first Should we calculate the correlation of the mean
#'     \code{design_perm} and \code{sv} (\code{calc_first = "mean"}), or
#'     should we calculate the mean of the correlations between
#'     \code{design_perm} and \code{sv} (\code{calc_first = "cor"})? This
#'     should only be changed by expert users.
#'
#' @export
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gale, David, and Lloyd S. Shapley. "College admissions and the stability of marriage." \emph{The American Mathematical Monthly} 69, no. 1 (1962): 9-15. \doi{10.1080/00029890.1962.11989827}.}
#'   \item{Gerard, D (2020). "Data-based RNA-seq simulations by binomial thinning." \emph{BMC Bioinformatics}. 21(1), 206. \doi{10.1186/s12859-020-3450-9}.}
#'   \item{Hornik K (2005). "A CLUE for CLUster Ensembles." \emph{Journal of Statistical Software}, 14(12). \doi{10.18637/jss.v014.i12}. \doi{10.18637/jss.v014.i12}.}
#'   \item{C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization: Algorithms and Complexity. Englewood Cliffs: Prentice Hall.}
#' }
#'
#' @return A matrix of correlations. The rows index the observed covariates
#'     and the columns index the surrogate variables. Element (i, j) is
#'     the estimated correlation between the ith variable in
#'     \code{design_perm} and the jth variable in \code{sv}.
#'
#' @examples
#' ## Generate the design matrices and set target correlation -----------------
#' n <- 10
#' design_perm <- cbind(rep(c(0, 1), each = n / 2),
#'                      rep(c(0, 1), length.out = n))
#' sv <- matrix(rnorm(n))
#' target_cor <- matrix(c(0.9, 0.1), ncol = 1)
#'
#' ## Get estimated true correlation ------------------------------------------
#' ## You should use a much larger iternum in practice
#' effective_cor(design_perm = design_perm,
#'               sv = sv,
#'               target_cor = target_cor,
#'               iternum = 10)
#'
effective_cor <- function(design_perm,
                          sv,
                          target_cor,
                          calc_first = c("cor", "mean"),
                          method = c("hungarian", "marriage"),
                          iternum = 1000) {
  ## Check input -------------------------------------------------------------
  calc_first <- match.arg(calc_first)
  assertthat::assert_that(is.matrix(design_perm))
  assertthat::assert_that(is.numeric(design_perm))
  assertthat::assert_that(is.matrix(sv))
  assertthat::assert_that(is.numeric(sv))
  assertthat::assert_that(is.matrix(target_cor))
  assertthat::assert_that(is.numeric(target_cor))
  assertthat::are_equal(nrow(sv), nrow(design_perm))
  assertthat::are_equal(ncol(sv), ncol(target_cor))
  assertthat::are_equal(ncol(design_perm), nrow(target_cor))
  assertthat::is.count(iternum)

  method <- match.arg(method)

  ## Get estimated correlation
  target_cor <- fix_cor(design_perm = design_perm, target_cor = target_cor)

  if (calc_first == "cor") {
    corarray <- array(NA, dim = c(ncol(design_perm), ncol(sv), iternum))
    ## latent_cor_array <- array(NA, dim = dim(corarray))
    for (index in seq_len(iternum)) {
      pout <- permute_design(design_perm = design_perm,
                             sv          = sv,
                             target_cor  = target_cor,
                             method      = method)
      corarray[, , index] <- stats::cor(pout$design_perm, sv)
      ## latent_cor_array[, , index] <- stats::cor(pout$latent_var, sv)
    }
    truecor <- apply(corarray, c(1, 2), mean)
  } else if (calc_first == "mean") {
    perm_array <- array(NA, dim = c(nrow(design_perm),
                                    ncol(design_perm),
                                    iternum))
    ## latent_array <- array(NA, dim = dim(perm_array))
    for (index in seq_len(iternum)) {
      pout <- permute_design(design_perm = design_perm,
                             sv          = sv,
                             target_cor  = target_cor,
                             method      = method)
      perm_array[, , index] <- pout$design_perm
      ## latent_array[, , index] <- pout$latent_var
    }
    truecor <- stats::cor(apply(perm_array, c(1, 2), mean), sv)
  }
  return(truecor)
}

#' Fixes an invalid target correlation.
#'
#' Shrinks the target correlation using a uniform scaling factor so that
#' the overall correlation matrix is positive semi-definite. The method
#' is described in detail in Gerard (2020).
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
#' @return A matrix of correlations the same dimension as \code{target_cor}.
#'    Actually, the returned matrix is \code{a * target_cor}, where \code{a}
#'    was determined to make the overall correlation matrix positive
#'    semi-definite.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D (2020). "Data-based RNA-seq simulations by binomial thinning." \emph{BMC Bioinformatics}. 21(1), 206. \doi{10.1186/s12859-020-3450-9}.}
#' }
#'
#' @export
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
  while (!valid) {
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
# @param mu A matrix of means. The rows index the independent samples,
#    the columns index the variables.
# @param sigma A covariance matrix of the columns.
rmvnorm <- function(mu, sigma) {
  stopifnot(nrow(sigma) == ncol(mu))
  cholout <- chol(sigma)
  simout <- matrix(stats::rnorm(n = prod(dim(mu))), nrow = nrow(mu)) %*% cholout + mu
  return(simout)
}
