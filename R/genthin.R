#########################
## General Thinning functions
#########################

#' Base binomial thinning function.
#'
#' Given a matrix of counts (\eqn{Y}) where \eqn{log_2(E[Y]) = Q},
#' a design matrix (\eqn{X}), and a matrix of coefficients (\eqn{B}),
#' \code{thin_diff} will generate a new matrix of counts such that
#' \eqn{log_2(E[Y]) = BX' + u1' + Q}, where \eqn{u} is some vector
#' of intercept coefficients. This function is used by all other
#' thinning functions. The method is
#' described in detail in Gerard (2019).
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
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_diff}}}{For the function most users should
#'       be using for general-purpose binomial thinning.}
#'   \item{\code{\link{thin_2group}}}{For the specific application of
#'       thinning in the two-group model.}
#'   \item{\code{\link{thin_lib}}}{For the specific application of
#'       library size thinning.}
#'   \item{\code{\link{thin_gene}}}{For the specific application of
#'       total gene expression thinning.}
#'   \item{\code{\link{thin_all}}}{For the specific application of
#'       thinning all counts uniformly.}
#' }
#'
#' @return A matrix of new RNA-seq read-counts. This matrix has the signal
#'     added from \code{designmat} and \code{coefmat}.
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
#' }
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
      stop(paste0("thin_base: tcrossprod(coefmat, designmat) produced positive entries\n",
                  "       and relative = FALSE. Either set relative = TRUE or change your\n",
                  "       coefficient and design matrices."))
    }
    qmat <- 2 ^ meanmat
  } else {
    qmat <- 2 ^ (meanmat - maxvec)
  }
  newmat      <- stats::rbinom(n = prod(dim(mat)), size = mat, prob = qmat)
  dim(newmat) <- dim(mat)
  return(newmat)
}

#' Binomial thinning for altering read-depth.
#'
#' Given a matrix of real RNA-seq counts, this function will apply a
#' thinning factor uniformly to every count in this matrix. This uniformly
#' lowers the read-depth for the entire dataset. The thinning factor should
#' be provided on the log2-scale. This is a specific application of the
#' binomial thinning approach in \code{\link{thin_diff}}. Though this particular
#' form of thinning was used by Robinson and Storey (2014) in the context
#' of deriving read-depth suggestions. It is also
#' described in detail in Gerard (2019).
#'
#' @inheritParams thin_diff
#' @param thinlog2 A numeric scalar. This is the amount to shrink each count
#'     in \code{mat} (on the log2-scale).  For
#'     example, a value of 0 means that we do not thin, a value of 1 means
#'     that we thin by a factor of 2, a value of 2 means we thin by a factor
#'     of 4, etc.
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#'   \item{\code{\link{thin_lib}}}{For thinning sample-wise.}
#'   \item{\code{\link{thin_gene}}}{For thinning gene-wise.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @inherit thin_diff return
#'
#' @author David Gerard
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
#'   \item{Robinson, David G., and John D. Storey. "subSeq: determining appropriate sequencing depth through efficient read subsampling." Bioinformatics 30, no. 23 (2014): 3424-3426.}
#' }
#'
#' @examples
#' ## Generate count data and set thinning factor
#' ## In practice, you would obtain mat from a real dataset, not simulate it.
#' set.seed(1)
#' n <- 10
#' p <- 1000
#' lambda <- 1000
#' mat <- matrix(lambda, ncol = n, nrow = p)
#' thinlog2 <- 1
#'
#' ## Thin read-depths
#' thout <- thin_all(mat = mat, thinlog2 = thinlog2)
#'
#' ## Compare empirical and theoretical proportions
#' mean(thout$mat) / lambda
#' 2 ^ -thinlog2
#'
thin_all <- function(mat, thinlog2) {
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  stopifnot(1 == length(thinlog2))
  assertthat::assert_that(thinlog2 > 0)

  thout <- thin_lib(mat      = mat,
                    thinlog2 = rep(thinlog2, ncol(mat)),
                    relative = FALSE)
  return(thout)
}

#' Binomial thinning for altering library size.
#'
#' Given a matrix of real RNA-seq counts, this function will apply a
#' separate, user-provided thinning factor to each sample. This uniformly
#' lowers the counts for all genes in a sample. The thinning factor
#' should be provided on the log2-scale. This is a specific application
#' of the binomial thinning approach in \code{\link{thin_diff}}. The method is
#' described in detail in Gerard (2019).
#'
#' @inheritParams thin_diff
#' @param thinlog2 A vector of numerics. Element i is the amount to thin
#'     (on the log2-scale) for sample i. For
#'     example, a value of 0 means that we do not thin, a value of 1 means
#'     that we thin by a factor of 2, a value of 2 means we thin by a factor
#'     of 4, etc.
#'
#' @inherit thin_diff return
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#'   \item{\code{\link{thin_gene}}}{For thinning gene-wise instead of
#'       sample-wise.}
#'   \item{\code{\link{thin_all}}}{For thinning all counts uniformly.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
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
#' lambda <- 1000
#' mat <- matrix(lambda, ncol = n, nrow = p)
#' thinlog2 <- rexp(n = n, rate = 1)
#'
#' ## Thin library sizes
#' thout <- thin_lib(mat = mat, thinlog2 = thinlog2)
#'
#' ## Compare empirical thinning proportions to specified thinning proportions
#' empirical_propvec <- colMeans(thout$mat) / lambda
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


#' Binomial thinning for altering total gene expression levels
#'
#' Given a matrix of real RNA-seq counts, this function will apply a
#' separate, user-provided thinning factor to each gene. This uniformly
#' lowers the counts for all samples in a gene. The thinning factor
#' should be provided on the log2-scale. This is a specific application
#' of the binomial thinning approach in \code{\link{thin_diff}}. The method is
#' described in detail in Gerard (2019).
#'
#'
#' @inheritParams thin_diff
#' @param thinlog2 A vector of numerics. Element i is the amount to thin
#'     (on the log2 scale) for gene i. For
#'     example, a value of 0 means that we do not thin, a value of 1 means
#'     that we thin by a factor of 2, a value of 2 means we thin by a factor
#'     of 4, etc.
#'
#' @inherit thin_diff return
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#'   \item{\code{\link{thin_lib}}}{For thinning sample-wise instead of
#'       gene-wise.}
#'   \item{\code{\link{thin_all}}}{For thinning all counts uniformly.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
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
#' lambda <- 1000
#' mat <- matrix(lambda, ncol = n, nrow = p)
#' thinlog2 <- rexp(n = p, rate = 1)
#'
#' ## Thin total gene expressions
#' thout <- thin_gene(mat = mat, thinlog2 = thinlog2)
#'
#' ## Compare empirical thinning proportions to specified thinning proportions
#' empirical_propvec <- rowMeans(thout$mat) / lambda
#' specified_propvec <- 2 ^ (-thinlog2)
#' plot(empirical_propvec, specified_propvec,
#'      xlab = "Empirical Thinning Proportion",
#'      ylab = "Specified Thinning Proportion")
#' abline(0, 1, col = 2, lwd = 2)
#'
thin_gene <- function(mat, thinlog2, relative = FALSE) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::are_equal(length(thinlog2), nrow(mat))
  thinlog2 <- c(thinlog2)
  assertthat::assert_that(is.numeric(thinlog2))
  stopifnot(thinlog2 >= 0)

  thout <- thin_diff(mat          = mat,
                     design_fixed = matrix(1, ncol = 1, nrow = ncol(mat)),
                     coef_fixed   = matrix(-thinlog2, ncol = 1),
                     relative     = relative)

  return(thout)

}

#' Binomial thinning in the two-group model.
#'
#' Given a matrix of real RNA-seq counts, this function will
#' randomly assign samples to one of two groups, draw
#' the log2-fold change in expression between two groups for each gene,
#' and add this signal to the RNA-seq counts matrix. The user may specify
#' the proportion of samples in each group, the proportion of null genes
#' (where the log2-fold change is 0),
#' and the signal function. This is a specific application of the
#' general binomial thinning approach implemented in \code{\link{thin_diff}}.
#'
#' The specific application of binomial thinning to the two-group model was
#' used in Gerard and Stephens (2017) and Gerard and Stephens (2018). This is
#' a specific case of the general method described in Gerard (2019).
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
#'     \code{i}th surrogate variable. The default is to set this to \code{NULL},
#'     in which case group assignment is made independently of any
#'     unobserved confounding.
#' @param alpha The scaling factor for the signal distribution from
#'     Stephens (2016). If \eqn{x_1, x_2, ..., x_n} are drawn from
#'     \code{signal_fun}, then the signal is set to
#'     \eqn{x_1 s_1^{\alpha}, x_2 s_2^{\alpha}, ..., x_n s_n^{\alpha}}, where
#'     \eqn{s_g} is the empirical standard deviation of gene \eqn{g}.
#'     Setting this to \code{0} means that the effects are exchangeable, setting
#'     this to \code{1} corresponds to the p-value prior of
#'     Wakefield (2009). You would rarely set this to anything but \code{0}
#'     or \code{1}.
#'
#' @inherit thin_diff return
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_diff}}}{For the more general thinning approach.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gale, David, and Lloyd S. Shapley. "College admissions and the stability of marriage." The American Mathematical Monthly 69, no. 1 (1962): 9-15.}
#'   \item{Gerard, David and Matthew Stephens (2017). "Unifying and generalizing methods for removing unwanted variation based on negative controls." \emph{arXiv} preprint arXiv:1705.08393.}
#'   \item{David Gerard and Matthew Stephens (2018). "Empirical Bayes shrinkage and false discovery rate estimation, allowing for unwanted variation." \emph{Biostatistics}, doi: \href{https://doi.org/10.1093/biostatistics/kxy029}{10.1093/biostatistics/kxy029}.}
#'   \item{Gerard David (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
#'   \item{Hansen, Ben B., and Stephanie Olsen Klopfer. "Optimal full matching and related designs via network flows." Journal of computational and Graphical Statistics 15, no. 3 (2006): 609-627.}
#'   \item{Hornik K (2005). "A CLUE for CLUster Ensembles." Journal of Statistical Software, 14(12). doi: 10.18637/jss.v014.i12}
#'   \item{C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization: Algorithms and Complexity. Englewood Cliffs: Prentice Hall.}
#'   \item{Stephens, Matthew. "False discovery rates: a new deal." Biostatistics 18, no. 2 (2016): 275-294.}
#'   \item{Wakefield, Jon. "Bayes factors for genome-wide association studies: comparison with P-values." Genetic epidemiology 33, no. 1 (2009): 79-86.}
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
                        corvec        = NULL,
                        alpha         = 0,
                        permute_method  = c("optmatch", "hungarian", "marriage")) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::are_equal(1L, length(prop_null))
  assertthat::are_equal(1L, length(signal_fun))
  assertthat::are_equal(1L, length(group_prop))
  assertthat::are_equal(1L, length(alpha))
  assertthat::assert_that(prop_null  <= 1,
                          prop_null  >= 0,
                          group_prop <= 1,
                          group_prop >= 0)
  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::assert_that(is.null(signal_params$n))
  assertthat::assert_that(is.numeric(alpha))

  ngene <- nrow(mat)
  nsamp <- ncol(mat)

  ## Generate coef ------------------------------------------------------------
  coef_perm <- matrix(0, nrow = ngene, ncol = 1)
  signal_params$n <- round(ngene * (1 - prop_null))
  if (signal_params$n > 0) {
    signal_vec <- c(do.call(what = signal_fun, args = signal_params))
    which_nonnull <- sort(sample(seq_len(ngene), size = signal_params$n))

    ## if alpha is non-zero
    if (abs(alpha) > 10 ^ -6) {
      matsd <- apply(X      = log2(mat[which_nonnull, ] + 0.5),
                     MARGIN = 1,
                     FUN    = stats::sd)
      signal_vec <- signal_vec * (matsd ^ alpha)
    }

    coef_perm[which_nonnull, 1] <- signal_vec
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
                     relative    = TRUE,
                     permute_method = permute_method)

  return(thout)
}

#' Binomial thinning for differential expression analysis.
#'
#' Given a matrix of real RNA-seq counts, this function will add a known
#' amount of signal to the count matrix. This signal is given in the form
#' of a Poisson / negative binomial / mixture of negative binomials
#' generalized linear model with a log (base 2) link. The user may
#' specify any arbitrary design matrix and coefficient matrix. The user
#' may also control for the amount of correlation between the observed
#' covariates and any unobserved surrogate variables. The method is
#' described in detail in Gerard (2019).
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
#'   \item{\eqn{M}}{Be a \eqn{G} by \eqn{N} of additional (unknown)
#'        unwanted variation.}
#' }
#' We assume that \eqn{Y} is Poisson distributed given \eqn{X_3} and
#' \eqn{Z} such that
#' \deqn{\log_2(EY) = \mu 1_N' + B_3X_3' + AZ' + M.}
#' \code{thin_diff()} will take as input \eqn{X_1}, \eqn{X_2}, \eqn{B_1},
#' \eqn{B_2}, and will output a \eqn{\tilde{Y}} and \eqn{W} such that
#' \eqn{\tilde{Y}} is Poisson distributed given \eqn{X_1}, \eqn{X_2}, \eqn{X_3},
#' \eqn{W}, \eqn{Z}, and \eqn{M} such that
#' \deqn{\log_2(E\tilde{Y}) \approx \tilde{\mu}1_N' + B_1X_1' + B_2X_2'W' + B_3X_3' + AZ' + M,}
#' where \eqn{W} is an \eqn{N} by \eqn{N} permutation matrix. \eqn{W} is randomly
#' drawn so that \eqn{WX_2} and \eqn{Z} are correlated approximately according
#' to the target correlation matrix.
#'
#' The Poisson assumption may be generalized to a mixture of negative binomials.
#'
#' @section Unestimable Components:
#'
#' It is possible to include an intercept term or a column from
#' \code{design_obs} into either \code{design_fixed} or \code{design_perm}.
#' This will not produce an error and the specified thinning will be applied.
#' However, If any column of \code{design_fixed} or
#' \code{design_perm} is a vector of ones or contains a column from
#' \code{design_obs}, then the corresponding columns in \code{coef_fixed}
#' or \code{coef_perm} cannot be estimated by \emph{any} method. This is
#' represented in the output by having duplicate columns in
#' \code{designmat} and \code{design_obs}.
#'
#' Including duplicate columns in \code{design_fixed} and \code{design_perm}
#' is also allowed but, again, will produce unestimable coefficients.
#'
#' Including an intercept term in \code{design_obs} will produce an error if
#' you are specifying correlated surrogate variables.
#'
#' @param mat A numeric matrix of RNA-seq counts. The rows index the genes and
#'     the columns index the samples.
#' @param design_fixed A numeric design matrix whose rows are fixed and not
#'     to be permuted. The rows index the samples and the columns index the
#'     variables. The intercept should \emph{not} be included
#'     (though see Section "Unestimable Components").
#' @param coef_fixed A numeric matrix. The coefficients corresponding to
#'     \code{design_fixed}. The rows index the genes and the columns index
#'     the variables.
#' @param design_perm A numeric design matrix whose rows are to be permuted
#'     (thus controlling the amount by which they are correlated with the
#'     surrogate variables). The rows index the samples and the columns index
#'     the variables. The intercept should \emph{not} be included
#'     (though see Section "Unestimable Components").
#' @param coef_perm A numeric matrix. The coefficients corresponding to
#'     \code{design_perm}. The rows index the genes and the columns index
#'     the variables.
#' @param target_cor A numeric matrix of target correlations between the
#'     variables in \code{design_perm} and the surrogate variables. The
#'     rows index the observed covariates and the columns index the surrogate
#'     variables. That is, \code{target_cor[i, j]} specifies the target
#'     correlation between the \code{i}th column of \code{design_perm} and the
#'     \code{j}th surrogate variable. The surrogate variables are estimated
#'     either using factor analysis or surrogate variable analysis (see the
#'     parameter \code{use_sva}).
#'     The number of columns in \code{target_cor} specifies the number of
#'     surrogate variables. Set \code{target_cor} to \code{NULL} to indicate
#'     that \code{design_perm} and the surrogate variables are independent.
#' @param use_sva A logical. Should we use surrogate variable analysis
#'     (Leek and Storey, 2008) using \code{design_obs}
#'     to estimate the hidden covariates (\code{TRUE})
#'     or should we just do an SVD on \code{log2(mat + 0.5)} after
#'     regressing out \code{design_obs} (\code{FALSE})? Setting this to
#'     \code{TRUE} allows the surrogate variables to be correlated with the
#'     observed covariates, while setting this to \code{FALSE} assumes that
#'     the surrogate variables are orthogonal to the observed covariates. This
#'     option only matters if \code{design_obs} is not \code{NULL}.
#'     Defaults to \code{FALSE}.
#' @param design_obs A numeric matrix of observed covariates that are NOT to
#'     be a part of the signal generating process. Only used in estimating the
#'     surrogate variables (if \code{target_cor} is not \code{NULL}).
#'     The intercept should \emph{not} be included (it will sometimes
#'     produce an error if it is included).
#' @param relative A logical. Should we apply relative thinning (\code{TRUE})
#'     or absolute thinning (\code{FALSE}). Only experts should change
#'     the default.
#' @param change_colnames A logical. Should we change the column-names
#'     of the design matrices (\code{TRUE}) or not (\code{FALSE})?
#'     Each new column name begins with either "O" (observed), "P" (permuted),
#'     or "F" (fixed), followed by a number. The letters correspond to
#'     whether the variables come from \code{design_obs}, \code{design_perm},
#'     or \code{design_fixed}. Setting this to \code{TRUE}
#'     also changes the column-names of the corresponding coefficient matrices.
#'     Defaults to \code{TRUE}.
#' @param permute_method Should we use the optimal matching technique from Hansen and
#'     Klopfer (2006) (\code{"optmatch"}), the Gale-Shapley algorithm
#'     for stable marriages (\code{"marriage"}) (Gale and Shapley, 1962)
#'     as implemented in the matchingR package, or the Hungarian algorithm
#'     (Papadimitriou and Steiglitz, 1982) (\code{"hungarian"})
#'     as implemented in the clue package (Hornik, 2005)?
#'     The \code{"optmatch"} method works really well
#'     but does take a lot more computational time if you have, say, 1000
#'     samples. If you use the \code{"optmatch"} option, you should note
#'     that the optmatch package uses a super strange license:
#'     \url{https://cran.r-project.org/package=optmatch/LICENSE}. If this
#'     license doesn't work for you (because you are not in academia, or
#'     because you don't believe in restrictive licenses), then
#'     try out the \code{"hungarian"} method.
#'
#' @return A list-like S3 object of class \code{ThinData}.
#' Components include some or all of the following:
#' \describe{
#'   \item{\code{mat}}{The modified matrix of counts.}
#'   \item{\code{designmat}}{The design matrix of variables used to simulate
#'       signal. This is made by column-binding \code{design_fixed} and the
#'       permuted version of \code{design_perm}.}
#'   \item{\code{coefmat}}{A matrix of coefficients corresponding to
#'       \code{designmat}.}
#'   \item{\code{design_obs}}{Additional variables that should be included in
#'       your design matrix in downstream fittings. This is made by
#'       column-binding the vector of 1's with \code{design_obs}.}
#'   \item{\code{sv}}{A matrix of estimated surrogate variables. In simulation
#'       studies you would probably leave this out and estimate your own
#'       surrogate variables.}
#'   \item{\code{cormat}}{A matrix of target correlations between the
#'       surrogate variables and the permuted variables in the design matrix.
#'       This might be different from the \code{target_cor} you input because
#'       we pass it through \code{\link{fix_cor}} to ensure
#'       positive semi-definiteness of the resulting covariance matrix.}
#'   \item{\code{matching_var}}{A matrix of simulated variables used to
#'       permute \code{design_perm} if the \code{target_cor} is not
#'       \code{NULL}.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{select_counts}}}{For subsampling the rows and columns
#'       of your real RNA-seq count matrix prior to applying binomial thinning.}
#'   \item{\code{\link{thin_2group}}}{For the specific application of
#'       \code{thin_diff} to the two-group model.}
#'   \item{\code{\link{thin_lib}}}{For the specific application of
#'       \code{thin_diff} to library size thinning.}
#'   \item{\code{\link{thin_gene}}}{For the specific application of
#'       \code{thin_diff} to total gene expression thinning.}
#'   \item{\code{\link{thin_all}}}{For the specific application of
#'       \code{thin_diff} to thinning all counts uniformly.}
#'   \item{\code{\link{thin_base}}}{For the underlying thinning function
#'       used in \code{thin_diff}.}
#'   \item{\code{\link[sva]{sva}}}{For the implementation of surrogate
#'       variable analysis.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gale, David, and Lloyd S. Shapley. "College admissions and the stability of marriage." The American Mathematical Monthly 69, no. 1 (1962): 9-15.}
#'   \item{Gerard D (2019). "Data-based RNA-seq Simulations by Binomial Thinning." \emph{bioRxiv}. doi: \href{https://doi.org/10.1101/758524}{10.1101/758524}.}
#'   \item{Hansen, Ben B., and Stephanie Olsen Klopfer. "Optimal full matching and related designs via network flows." Journal of computational and Graphical Statistics 15, no. 3 (2006): 609-627.}
#'   \item{Hornik K (2005). "A CLUE for CLUster Ensembles." Journal of Statistical Software, 14(12). doi: 10.18637/jss.v014.i12}
#'   \item{Leek, Jeffrey T., and John D. Storey. "A general framework for multiple testing dependence." Proceedings of the National Academy of Sciences 105, no. 48 (2008): 18718-18723.}
#'   \item{C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization: Algorithms and Complexity. Englewood Cliffs: Prentice Hall.}
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
#' plot(thout$coefmat[, 1], coef_est[, 1],
#'      main = "First Variable",
#'      xlab = "Coefficient",
#'      ylab = "Estimated Coefficient")
#' abline(0, 1, col = 2, lwd = 2)
#'
#' plot(thout$coefmat[, 2], coef_est[, 2],
#'      main = "Second Variable",
#'      xlab = "Coefficient",
#'      ylab = "Estimated Coefficient")
#' abline(0, 1, col = 2, lwd = 2)
#'
#' ## But estimated coefficient of the first variable is slightly closer when
#' ## the surrogate variable is included.
#' Ynew <- log2(t(thout$mat) + 0.5)
#' X <- cbind(thout$designmat, thout$sv)
#' coef_est <- t(coef(lm(Ynew ~ X))[2:3, ])
#'
#' plot(thout$coefmat[, 1], coef_est[, 1],
#'      main = "First Variable",
#'      xlab = "Coefficient",
#'      ylab = "Estimated Coefficient")
#' abline(0, 1, col = 2, lwd = 2)
#'
#' plot(thout$coefmat[, 2], coef_est[, 2],
#'      main = "Second Variable",
#'      xlab = "Coefficient",
#'      ylab = "Estimated Coefficient")
#' abline(0, 1, col = 2, lwd = 2)
#'
thin_diff <- function(mat,
                      design_fixed    = NULL,
                      coef_fixed      = NULL,
                      design_perm     = NULL,
                      coef_perm       = NULL,
                      target_cor      = NULL,
                      use_sva         = FALSE,
                      design_obs      = NULL,
                      relative        = TRUE,
                      change_colnames = TRUE,
                      permute_method  = c("optmatch", "hungarian", "marriage")) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::assert_that(is.logical(use_sva))
  assertthat::assert_that(is.logical(relative))
  assertthat::assert_that(is.logical(change_colnames))
  assertthat::are_equal(1L, length(use_sva))
  assertthat::are_equal(1L, length(relative))
  assertthat::are_equal(1L, length(change_colnames))
  stopifnot(mat >= 0)

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

  assertthat::assert_that(is.matrix(design_fixed))
  assertthat::assert_that(is.matrix(design_perm))
  assertthat::assert_that(is.matrix(design_obs))
  assertthat::assert_that(is.matrix(coef_fixed))
  assertthat::assert_that(is.matrix(coef_perm))
  assertthat::assert_that(is.numeric(design_fixed))
  assertthat::assert_that(is.numeric(design_perm))
  assertthat::assert_that(is.numeric(design_perm))
  assertthat::assert_that(is.numeric(coef_fixed))
  assertthat::assert_that(is.numeric(coef_perm))
  assertthat::are_equal(ncol(mat), nrow(design_fixed), nrow(design_perm), nrow(design_obs))
  assertthat::are_equal(nrow(mat), nrow(coef_fixed), nrow(coef_perm))
  assertthat::are_equal(ncol(design_fixed), ncol(coef_fixed))
  assertthat::are_equal(ncol(design_perm), ncol(coef_perm))

  permute_method <- match.arg(permute_method)
  if (!is.null(target_cor) & permute_method == "optmatch" & !requireNamespace("optmatch", quietly = TRUE)) {
    stop(paste0("\nPackage optmatch must be installed to use `permute_method = \"optmatch\"`\n",
                "You can install it with\n\n",
                "install.packages(\"optmatch\")\n\n",
                "Note that optmatch uses a strange non-standard license:\n",
                "https://cran.r-project.org/package=optmatch/LICENSE\n"))
  } else if (!is.null(target_cor) & permute_method == "optmatch") {
    message_fun("optmatch")
  }

  ## Permute ------------------------------------------------------------------
  if (is.null(target_cor) | ncol(design_perm) == 0) {
    design_perm <- design_perm[sample(seq_len(nrow(design_perm))), , drop = FALSE]
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
                 design_obs   = design_obs,
                 use_sva      = use_sva)

    ## Permute design matrix ----------------------------
    pout <- permute_design(design_perm = design_perm,
                           sv          = sv,
                           target_cor  = new_cor,
                           method      = permute_method)
    design_perm <- pout$design_perm
    latent_var  <- pout$latent_var
  }

  ## Fix column names ---------------------------------------------------------
  if (change_colnames) {
    if (ncol(design_fixed) > 0) {
      colnames(design_fixed) <- paste0("F", seq_len(ncol(design_fixed)))
      colnames(coef_fixed)   <- paste0("F", seq_len(ncol(design_fixed)))
    }
    if (ncol(design_perm) > 0) {
      colnames(design_perm)   <- paste0("P", seq_len(ncol(design_perm)))
      colnames(coef_perm)     <- paste0("P", seq_len(ncol(design_perm)))
    }
    if (ncol(design_obs) > 0) {
      colnames(design_obs)   <- paste0("O", seq_len(ncol(design_obs)))
    }
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
                 designmat    = cbind(designmat),
                 coefmat      = coefmat,
                 design_obs   = cbind("(Intercept)" = 1, design_obs),
                 sv           = sv,
                 cormat       = new_cor,
                 matching_var = latent_var)

  class(retval) <- "ThinData"

  return(retval)
}
