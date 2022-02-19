## Generate counts from the theoretical model.

#' Simulate RNA-seq data from the (unrealistic) theoretical model.
#'
#' The main difference between \code{seqgen_base} and \code{seqgen_diff}
#' is the output format.
#'
#' You should be using \code{\link{thin_base}} instead of this. If you must
#' simulate from the theoretical model, I would recommend trying out
#' the powsimR package (\url{https://github.com/bvieth/powsimR}) from
#' Vieth et al. (2017).
#'
#' @inheritParams seqgen_diff
#'
#' @return A matrix of RNA-seq counts. The rows index the genes and the columns
#'     index the samples.
#'
#' @references
#' \itemize{
#'   \item{Vieth, Beate, Christoph Ziegenhain, Swati Parekh, Wolfgang Enard, and Ines Hellmann. "powsimR: power analysis for bulk and single cell RNA-seq experiments." \emph{Bioinformatics} 33, no. 21 (2017): 3486-3488. \doi{10.1093/bioinformatics/btx435}.}
#' }
#'
#' @author David Gerard
#' @noRd
seqgen_base <- function(designmat,
                        coefmat,
                        dispvec) {

  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(designmat))
  assertthat::assert_that(is.matrix(coefmat))
  assertthat::assert_that(is.numeric(designmat))
  assertthat::assert_that(is.numeric(coefmat))
  assertthat::assert_that(is.numeric(dispvec))
  assertthat::assert_that(all(dispvec > 0))
  assertthat::are_equal(nrow(coefmat), length(dispvec))
  assertthat::are_equal(ncol(coefmat), ncol(designmat))

  ngene <- nrow(coefmat)
  nsamp <- nrow(designmat)

  ## Get mean -----------------------------------------------------------------
  mumat <- 2 ^ tcrossprod(coefmat, designmat)

  ## Get draws ----------------------------------------------------------------
  mat <- stats::rnbinom(n    = ngene * nsamp,
                        size = 1 / dispvec,
                        mu   = mumat)
  dim(mat) <- dim(mumat)

  return(mat)
}

#' Theoretical sequence simulation for differential expression analysis.
#'
#' Generate a matrix of RNA-seq counts from a theoretical negative-binomial
#' model.
#'
#' You should be using \code{\link{thin_diff}} instead of this. If you must
#' simulate from the theoretical model, I would recommend trying out
#' the powsimR package (\url{https://github.com/bvieth/powsimR}) from
#' Vieth et al. (2017).
#'
#' @param designmat A numeric design matrix. The rows index the samples
#'     and the columns index the variables.
#' @param coefmat A numeric matrix of coefficients. The rows index the
#'     genes and the columns index the variables.
#' @param dispvec A vector of negative binomial dispersion parameters.
#'     Element i is the dispersion for gene i.
#' @param design_sv A numeric design matrix for the surrogate variables.
#'     The rows index the samples and the columns index the surrogate
#'     variables.
#' @param coef_sv A numeric matrix of coefficients for the surrogate
#'     variables. The rows index the genes and the columns index the
#'     variables.
#'
#' @references
#' \itemize{
#'   \item{Vieth, Beate, Christoph Ziegenhain, Swati Parekh, Wolfgang Enard, and Ines Hellmann. "powsimR: power analysis for bulk and single cell RNA-seq experiments." \emph{Bioinformatics} 33, no. 21 (2017): 3486-3488. \doi{10.1093/bioinformatics/btx435}.}
#' }
#'
#' @return A list object with some or all of the following elements:
#' \describe{
#'  \item{\code{mat}}{A matrix of RNA-seq counts. The rows index the genes and
#'      the columns index the samples.}
#'  \item{\code{designmat}}{The user-provided \code{designmat}.}
#'  \item{\code{coefmat}}{The user-provided \code{coefmat}.}
#'  \item{\code{sv}}{The user-provided \code{design_sv}.}
#'  \item{\code{coef_sv}}{The user-provided \code{coef_sv}.}
#' }
#'
#' @author David Gerard
#' @noRd
seqgen_diff <- function(designmat,
                        coefmat,
                        dispvec,
                        design_sv = NULL,
                        coef_sv = NULL) {
  ## Check input --------------------------------------------------------------
  ngene <- nrow(coefmat)
  nsamp <- nrow(designmat)
  if (is.null(design_sv) | is.null(coef_sv)) {
    design_sv <- matrix(ncol = 0L, nrow = nsamp)
    class(design_sv) <- "numeric"
    coef_sv <- matrix(ncol = 0L, nrow = ngene)
    class(coef_sv) <- "numeric"
  }
  assertthat::assert_that(is.matrix(designmat))
  assertthat::assert_that(is.matrix(coefmat))
  assertthat::assert_that(is.matrix(design_sv))
  assertthat::assert_that(is.matrix(coef_sv))
  assertthat::assert_that(is.numeric(designmat))
  assertthat::assert_that(is.numeric(coefmat))
  assertthat::assert_that(is.numeric(dispvec))
  assertthat::assert_that(is.numeric(design_sv))
  assertthat::assert_that(is.numeric(coef_sv))
  assertthat::are_equal(nrow(coefmat), nrow(coef_sv), length(dispvec))
  assertthat::are_equal(nrow(designmat), nrow(design_sv))
  assertthat::are_equal(ncol(designmat), ncol(coefmat))
  assertthat::are_equal(ncol(design_sv), ncol(coef_sv))
  assertthat::assert_that(all(dispvec > 0))

  ## Generate RNA-seq ---------------------------------------------------------
  mat <- seqgen_base(designmat = cbind(designmat, design_sv),
                     coefmat   = cbind(coefmat, coef_sv),
                     dispvec   = dispvec)

  ## Return -------------------------------------------------------------------
  retlist <- list(mat       = mat,
                  designmat = designmat,
                  coefmat   = coefmat,
                  sv        = design_sv,
                  coef_sv   = coef_sv)

  return(retlist)
}


#' Theoretical sequence simulation in the two-group model.
#'
#' Generate a matrix of RNA-seq counts from a theoretical negative-binomial
#' model. This is in the special case of the two-group model.
#'
#' You should be using \code{\link{thin_2group}} instead of this. If you must
#' simulate from the theoretical model, I would recommend trying out
#' the powsimR package (\url{https://github.com/bvieth/powsimR}) from
#' Vieth et al. (2017).
#'
#' @param ngene The number of genes.
#' @param nsamp The number of samples.
#' @param dispvec A numeric vector of dispersions.
#' @param prop_null The proportion of genes that are null.
#' @param signal_fun The signal function for the non-null genes.
#'     A function that takes at least the argument \code{n}
#'     and returns \code{n} draws. Additional parameters may be passed
#'     through \code{signal_params}.
#' @param signal_params Additional arguments to be passed to \code{signal_fun}.
#' @param intercept_fun The function for the intercept coefficients. A function
#'     that takes at least the argument \code{n} and returns \code{n} draws.
#'     Additional parameters may be passed through \code{intercept_params}.
#' @param intercept_params Additional arguments to be passed to
#'     \code{intercept_fun}.
#' @param libvec A vector of library size multiplicative factors (not
#'     on the log2 scale).
#' @param group_prop The proportion of samples to be placed in group 1.
#' @param design_sv Optional numeric design matrix of surrogate variables. The
#'     rows index the samples and the columns index the surrogate variables.
#' @param coef_sv Optional numeric matrix of coefficients of the surrogate
#'     variables. The rows index the genes and the columns index the surrogate
#'     variables.
#'
#' @references
#' \itemize{
#'   \item{Vieth, Beate, Christoph Ziegenhain, Swati Parekh, Wolfgang Enard, and Ines Hellmann. "powsimR: power analysis for bulk and single cell RNA-seq experiments." \emph{Bioinformatics} 33, no. 21 (2017): 3486-3488. \doi{10.1093/bioinformatics/btx435}.}
#' }
#'
#' @seealso
#' \describe{
#'   \item{\code{\link[DESeq2]{makeExampleDESeqDataSet}}}{For a very similar
#'       function from the DESeq2 package.}
#'   \item{\code{\link{seqgen_diff}}}{For the underlying sequence simulation
#'       function.}
#' }
#'
#' @inherit seqgen_diff return
#'
#' @author David Gerard
#' @noRd
seqgen_2group  <- function(ngene,
                           nsamp,
                           dispvec          = rep(0.1, times = ngene),
                           prop_null        = 1,
                           signal_fun       = stats::rnorm,
                           signal_params    = list(mean = 0, sd = 1),
                           intercept_fun    = stats::rnorm,
                           intercept_params = list(mean = 4, sd = 2),
                           libvec           = rep(1, times = nsamp),
                           group_prop       = 0.5,
                           design_sv        = NULL,
                           coef_sv          = NULL) {
  ## Check input --------------------------------------------------------------
  assertthat::is.count(ngene)
  assertthat::is.count(nsamp)
  if (is.null(design_sv) | is.null(coef_sv)) {
    design_sv <- matrix(ncol = 0L, nrow = nsamp)
    class(design_sv) <- "numeric"
    coef_sv <- matrix(ncol = 0L, nrow = ngene)
    class(coef_sv) <- "numeric"
  }
  assertthat::assert_that(is.numeric(dispvec))
  assertthat::assert_that(is.numeric(prop_null))
  assertthat::assert_that(is.numeric(group_prop))
  assertthat::assert_that(is.numeric(design_sv))
  assertthat::assert_that(is.numeric(coef_sv))
  assertthat::assert_that(is.matrix(design_sv))
  assertthat::assert_that(is.matrix(coef_sv))
  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.function(intercept_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::assert_that(is.list(intercept_params))
  assertthat::are_equal(1, length(prop_null), length(group_prop))
  assertthat::assert_that(prop_null >= 0, prop_null <= 1)
  assertthat::assert_that(group_prop >= 0, group_prop <= 1)
  assertthat::are_equal(ngene, length(dispvec), nrow(coef_sv))
  assertthat::are_equal(nsamp, length(libvec), nrow(design_sv))
  assertthat::assert_that(all(dispvec > 0))
  assertthat::assert_that(all(libvec >= 0))
  assertthat::assert_that(is.null(signal_params$n))
  assertthat::assert_that(is.null(intercept_params$n))

  ## Generate designmat and coefmat -------------------------------------------
  n1 <- round(group_prop * nsamp)
  groupvec <- rep(0, times = nsamp)
  groupvec[sample(seq_along(groupvec), size = n1)] <- 1
  designmat <- cbind("(Intercept)" = 1, "Treatmemt" = groupvec)

  intercept_params$n <- ngene
  beta0vec <- do.call(what = intercept_fun, args = intercept_params)

  signal_params$n <- ngene
  beta1vec <- do.call(what = signal_fun, args = signal_params)

  coefmat <- cbind("(Intercept)" = beta0vec, "Treatment" = beta1vec)

  ## Simulate -----------------------------------------------------------------
  retlist <- seqgen_diff(designmat = designmat,
                         coefmat   = coefmat,
                         dispvec   = dispvec,
                         design_sv = cbind(design_sv, "Library" = log2(libvec)),
                         coef_sv   = cbind(coef_sv, "Library" = 1))

  return(retlist)
}
