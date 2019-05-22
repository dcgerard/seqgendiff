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
#' ## In practice, you would obtain Y from a real dataset, not simulate it.
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
#' plot(B, Bhat, xlab = "Coefficients", ylab = "Coefficient Estimates")
#' abline(0, 1, col = 2, lwd = 2)
#'
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
#' Given a matrix of real RNA-seq counts, this function will apply a
#' separate, user-provided thinning factor to each sample. This uniformly
#' lowers the counts for all genes in a sample. The thinning factor
#' should be provided on the log2-scale. This is a specific application
#' of the Poisson thinning approach in \code{\link{thin_diff}}.
#'
#' @inheritParams thin_diff
#' @param thinlog2 A vector of numerics. Element i is the amount to thin (on the log2 scale). For
#'     example, a value of 0 means that we do not thin, a value of 1 means
#'     that we thin by a factor of 2, a value of 2 means we thin by a factor
#'     of 4, etc.
#'
#' @inherit thin_diff return
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @examples
#' ## Generate count data and thinning factors
#' ## In practice, you would obtain mat from a real dataset, not simulate it.
#' set.seed(1)
#' n <- 10
#' p <- 1000
#' mat <- matrix(1000, ncol = n, nrow = p)
#' thinlog2 <- rexp(n = n, rate = 1)
#'
#' ## Thin library sizes
#' thout <- thin_lib(mat = mat, thinlog2 = thinlog2)
#'
#' ## Compare empirical thinning proportions to specified thinning proportions
#' empirical_propvec <- colMeans(thout$mat) / 1000
#' specified_propvec <- 2 ^ (-thinlog2)
#' empirical_propvec
#' specified_propvec
#'
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
#' Given a matrix of real RNA-seq counts, this function will
#' randomly assign samples to one of two groups, draw
#' the log2-fold change in expression between two groups for each gene,
#' and add this signal to the RNA-seq counts matrix. The user may specify
#' the proportion of samples in each group, the proportion of null genes
#' (where the log2-fold change is 0),
#' and the signal function. This is a specific application of the
#' general Poisson thinning approach implemented in \code{\link{thin_diff}}.
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
#' @inherit thin_diff return
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#' }
#'
#' @examples
#' ## Simulate data from given matrix of counts
#' ## In practice, you would obtain Y from a real dataset, not simulate it.
#' set.seed(1)
#' nsamp <- 10
#' ngene <- 1000
#' Y <- matrix(stats::rpois(nsamp * ngene, lambda = 50), nrow = ngene)
#' thinout <- thin_2group(mat           = Y,
#'                        prop_null     = 0.9,
#'                        signal_fun    = stats::rexp,
#'                        signal_params = list(rate = 0.5))
#'
#' ## 90 percent of genes are null
#' mean(abs(thinout$coef) < 10^-6)
#'
#' ## Check the estimates of the log2-fold change
#' Ynew <- log2(t(thinout$mat + 0.5))
#' X    <- thinout$designmat
#' Bhat <- coef(lm(Ynew ~ X))["X", ]
#' plot(thinout$coefmat,
#'      Bhat,
#'      xlab = "log2-fold change",
#'      ylab = "Estimated log2-fold change")
#' abline(0, 1, col = 2, lwd = 2)
#'
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
                     target_cor  = target_cor,
                     relative    = TRUE)

  return(thout)
}

#' Poisson thinning for differential expression analysis.
#'
#' Given a matrix of real RNA-seq counts, this function will add a known
#' amount of signal to the count matrix. This signal is given in the form
#' of a Poisson generalized linear model with a log (base 2) link. The user may
#' specify any arbitrary design matrix and coefficient matrix. The user
#' may also control for the amount of correlation between the observed
#' covariates and any unobserved surrogate variables.
#'
#'
#' @section Mathematical Formulation:
#' Let
#' \describe{
#'   \item{\eqn{N}}{Be the number of samples.}
#'   \item{\eqn{G}}{Be the number of genes.}
#'   \item{\eqn{Y}}{Be an \eqn{G} by \eqn{N} matrix of real RNA-seq counts.
#'       This is \code{mat}.}
#'   \item{\eqn{X_1}}{Be an \eqn{N} by \eqn{P_1} user-provided design matrix.
#'       This is \code{design_fixed}.}
#'   \item{\eqn{X_2}}{Be an \eqn{N} by \eqn{P_2} user-provided design matrix.
#'       This is \code{design_perm}.}
#'   \item{\eqn{X_3}}{Be an \eqn{N} by \eqn{P_3} matrix of known covariates.
#'       This is \code{design_obs}.}
#'   \item{\eqn{Z}}{Be an \eqn{N} by \eqn{K} matrix of unobserved surrogate
#'        variables. This is estimated when \code{target_cor} is not
#'        \code{NULL}.}
#' }
#' We assume that \eqn{Y} is Poisson distributed given \eqn{X_3} and
#' \eqn{Z} such that
#' \deqn{\log_2(EY) = 1_N\mu' + B_3X_3' + AZ'.}
#' \code{thin_diff()} will take as inpute \eqn{X_1}, \eqn{X_2}, \eqn{B_1},
#' \eqn{B_2}, and will output a \eqn{\tilde{Y}} and \eqn{W} such that
#' \eqn{\tilde{Y}} is Poisson distributed given \eqn{X_1}, \eqn{X_2}, \eqn{X_3},
#' \eqn{W}, and \eqn{Z} such that
#' \deqn{\log_2(E\tilde{Y}) \approx 1_N\tilde{\mu}' + B_1X_1' + B_2X_2'W' + B_3X_3' + AZ',}
#' where \eqn{W} is an \eqn{N} by \eqn{N} permutation matrix. \eqn{W} is randomly
#' drawn so that \eqn{X_2} and \eqn{Z} are correlated approximately according
#' to a target correlation matrix.
#'
#'
#' @param mat A numeric matrix of counts. The rows index the genes and the
#'     columns index the samples (as is usual in RNA-seq).
#' @param design_fixed A numeric design matrix whose rows are fixed and not permuted.
#'     The rows index the samples and the columns index the variables.
#'     The intercept should \emph{not} be included.
#' @param coef_fixed A numeric matrix. The coefficients corresponding to \code{design_fixed}.
#'     The rows index the genes and the columns index the variables.
#' @param design_perm A numeric design matrix whose rows are to be permuted (thus
#'     controlling the amount by which they are correlated with the surrogate variables).
#'     The rows index the samples and the columns index the variables.
#'     The intercept should \emph{not} be included.
#' @param coef_perm A numeric matrix. The coefficients corresponding to \code{design_perm}.
#'     The rows index the genes and the columns index the variables.
#' @param target_cor A numeric matrix of target correlations betweens the variables in
#'     \code{design_perm} and the surrogate variables. The rows index the
#'     observed covariates and the columns index the surrogate variables.
#'     The number of columns indicates the number of surrogate variables.
#'     Set this to \code{NULL} to indicate uncorrelated.
#' @param use_sva A logical. Should we use SVA using \code{design_fixed}
#'     and \code{design_obs} to estimate the hidden covariates (\code{TRUE})
#'     or should we just do an SVD on \code{log2(mat + 0.5)} after
#'     subtracing row means (\code{FALSE})? Defaults to \code{FALSE}.
#' @param design_obs A numeric matrix of observed covariates that we are NOT to be a
#'     part of the signal generating process. Only used in estimating the
#'     surrogate variables if \code{use_sva = TRUE}.
#'     The intercept should \emph{not} be included.
#' @param relative A logical. Should we apply relative thinning (\code{TRUE})
#'     or absolute thinning (\code{FALSE}). Only experts should change
#'     the default.
#'
#' @return An S3 object of class \code{ThinData}. Components include some or
#' all of the following:
#' \describe{
#'   \item{\code{mat}}{The modified matrix of counts.}
#'   \item{\code{designmat}}{The design matrix, excluding an intercept term.}
#'   \item{\code{coefmat}}{A matrix of coefficients corresponding to \code{designmat}}
#'   \item{\code{sv}}{A matrix of estimated surrogate variables.}
#'   \item{\code{cormat}}{A matrix of target correlations between the
#'       surrogate variables and the permuted variables in the design matrix.}
#'   \item{\code{matching_var}}{A matrix of simulated variables used to
#'       permute the permuted components of the design matrix.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{thin_2group}}}{For the specific application of
#'       \code{thin_diff} to the two-group model.}
#'   \item{\code{\link{thin_lib}}}{For the specific application of
#'       \code{thin_diff} to library size thinning.}
#'   \item{\code{\link{thin_base}}}{For the underlying thinning function
#'       used in \code{thin_diff}.}
#' }
#'
#'
#' @examples
#' ## Generate simulated data with surrogate variables
#' ## In practice, you would obtain mat from a real dataset, not simulate it.
#' set.seed(1)
#' n <- 10
#' p <- 1000
#' Z <- matrix(abs(rnorm(n, sd = 4)))
#' alpha <- matrix(abs(rnorm(p, sd = 1)))
#' mat <- round(2^(alpha %*% t(Z) + abs(matrix(rnorm(n * p, sd = 5),
#'                                             nrow = p,
#'                                             ncol = n))))
#'
#' ## Choose simulation parameters
#' design_perm <- cbind(rep(c(0, 1), length.out = n), runif(n))
#' coef_perm <- matrix(rnorm(p * ncol(design_perm), sd = 6), nrow = p)
#'
#' ## Specify one surrogate variable (number of columns in taget_cor),
#' ## highly correlated with first observed covariate and uncorrelated
#' ## with second observed covariate
#' target_cor <- matrix(c(0.9, 0))
#'
#' ## Thin
#' thout <- thin_diff(mat = mat,
#'                    design_perm = design_perm,
#'                    coef_perm = coef_perm,
#'                    target_cor = target_cor)
#'
#' ## target_cor approximates correlation between estimated surrogate variable
#' ## and matching variable.
#' cor(thout$matching_var, thout$sv)
#'
#' ## Estimated surrogate variable is associated with true surrogate variable
#' ## (because the signal is strong in this case)
#' plot(Z, thout$sv, xlab = "True SV", ylab = "Estimated SV")
#'
#' ## So target_cor approximates correlation between surrogate variable and
#' ## matching variables
#' cor(thout$matching_var, Z)
#'
#' ## Correlation between permuted covariates and surrogate variables are less
#' ## close to target_cor
#' cor(thout$designmat, Z)
#'
#' ## Estimated signal is correlated to true single. First variable is slightly
#' ## biased because the surrogate variable is not included.
#' Ynew <- log2(t(thout$mat) + 0.5)
#' X <- thout$designmat
#' coef_est <- t(coef(lm(Ynew ~ X))[2:3, ])
#'
#' plot(thout$coefmat[, 1], coef_est[, 1])
#' abline(0, 1, col = 2, lwd = 2)
#'
#' plot(thout$coefmat[, 2], coef_est[, 2])
#' abline(0, 1, col = 2, lwd = 2)
#'
#' ## But estimated coefficient of the first variable is slightly closer when
#' ## the surrogate variable is included.
#' Ynew <- log2(t(thout$mat) + 0.5)
#' X <- cbind(thout$designmat, thout$sv)
#' coef_est <- t(coef(lm(Ynew ~ X))[2:3, ])
#'
#' plot(thout$coefmat[, 1], coef_est[, 1])
#' abline(0, 1, col = 2, lwd = 2)
#'
#' plot(thout$coefmat[, 2], coef_est[, 2])
#' abline(0, 1, col = 2, lwd = 2)
#'
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

    ## Estimate surrogate variables ---------------------
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
                 designmat    = cbind(designmat, design_obs),
                 coefmat      = coefmat,
                 sv           = sv,
                 cormat       = new_cor,
                 matching_var = latent_var)

  class(retval) <- "ThinData"

  return(retval)
}


