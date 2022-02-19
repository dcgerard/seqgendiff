
#' Apply Poisson thinning to a matrix of count data.
#'
#' This is now defunct. Please try out \code{\link{select_counts}} and
#' \code{\link{thin_2group}}.
#'
#' Given a matrix of RNA-seq counts, this function will randomly select two groups of
#' samples and add signal to a known proportion of the genes. This signal
#' is the log (base 2) effect size of the group indicator in a linear model.
#' The user may specify the distribution of the effects.
#'
#' The Poisson thinning approach first randomly assigns samples to be in one of two groups. Then,
#' given this assignment, will Binomially sample counts with a sample size of the gene expression
#' counts and a probability that is a function of the effect size. For details, see
#' Gerard and Stephens (2021).
#'
#' @param mat A matrix of count data. The rows index the individuals and
#'     the columns index the genes.
#' @param nsamp The number of samples to select from \code{mat}.
#' @param ngene The number of genes to select from \code{mat}.
#' @param gselect How should we select the subset of genes? Should we choose
#'     the \code{ngene} most median expressed genes (\code{"max"}), a random sample
#'     of the genes (\code{"random"}), a random sample of the most expressed
#'     genes (\code{"rand_max"}), a user-provided list (\code{"custom"}), or by maximum
#'     mean expression level (\code{"mean_max"})?
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
#'     Not used if \code{gselect = "custom"}.
#' @param prop_null The proportion of genes that are null.
#' @param alpha If \eqn{b} is an effect and \eqn{s} is an empirical standard deviation, then
#'     we model \eqn{b/s^\alpha} as being exchangeable.
#' @param group_assign How should we assign groups? Exactly specifying the
#'     proportion of individuals in each group (\code{"frac"}), with a
#'     Bernoulli distribution (\code{"random"}), or correlated with latent factors
#'     (\code{"cor"})? If \code{group_assign = "cor"}, then you have to specify
#'     \code{corvec}. If \code{group_assign = "frac"} or
#'     \code{group_assign = "random"}, then the proportion of samples in each
#'     group is specified with the \code{group_prop} argument.
#' @param group_prop The proportion of individuals that are in group 1.
#'     This proportion is deterministic if \code{group_assign = "frac"}, and
#'     is the expected proportion if \code{group_assign = "random"}. This
#'     argument is not used if \code{group_assign = "cor"}.
#' @param corvec A vector of correlations. \code{corvec[i]} is the correlation
#'     of the latent group assignment vector with the ith latent confounder.
#'     Only used if \code{group_assign = "cor"}. This vector is constrained
#'     so that \code{crossprod(corvec) < 1}. The number of latent factors
#'     is taken to be the length of corvec. Note that the correlations of the
#'     latent factors with the observed group-assignment vector (instead of the
#'     latent group-assignment vector) will be \code{corvec * sqrt(2 / pi)}.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item{\code{Y}: }{A matrix of altered counts with \code{nsamp} rows
#'        and \code{ngene} columns.}
#'  \item{\code{X}: }{A design matrix. The first column contains a vector ones (for an
#'        intercept term) and the second column contains an indicator for group membership.}
#'  \item{\code{beta}: }{The approximately true effect sizes of \eqn{log(Y) ~ X\beta}.}
#'  \item{\code{corassign}: }{The output from the call to \code{\link{corassign}}.
#'        Only returned if \code{group_assign = "cor"}.}
#' }
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., and Stephens, M. (2021). "Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls." \emph{Statistica Sinica}, 31(3), 1145-1166 \doi{10.5705/ss.202018.0345}.}
#' }
#'
#' @examples
#' ## Simulate data from given matrix of counts
#' ## In practice, you would obtain Y from a real dataset, not simulate it.
#' set.seed(1)
#' nsamp <- 10
#' ngene <- 1000
#' Y <- matrix(stats::rpois(nsamp * ngene, lambda = 50), nrow = ngene)
#'
#'
#' ## Apply thinning
#' poisout <- poisthin(mat           = t(Y),
#'                     nsamp         = 9,
#'                     ngene         = 999,
#'                     signal_fun    = stats::rnorm,
#'                     signal_params = list(mean = 0, sd = 1),
#'                     prop_null     = 0.9)
#'
#' ## Dimension of count matrix is smaller.
#' dim(poisout$Y)
#'
#' ## Can verify signal was added by estimating it with lm().
#' betahat <- coef(lm(log2(poisout$Y + 1) ~ poisout$X[, 2]))[2, ]
#' plot(poisout$beta, betahat, xlab = "Coefficients", ylab = "Estimates")
#' abline(0, 1, col = 2, lty = 2)
#'
#' @export
poisthin <- function(mat,
                     nsamp         = nrow(mat),
                     ngene         = ncol(mat),
                     gselect       = c("max",
                                       "random",
                                       "rand_max",
                                       "custom",
                                       "mean_max"),
                     gvec          = NULL,
                     skip_gene     = 0L,
                     signal_fun    = stats::rnorm,
                     signal_params = list(mean = 0, sd = 1),
                     prop_null     = 1,
                     alpha         = 0,
                     group_assign  = c("frac",
                                       "random",
                                       "cor"),
                     group_prop    = 0.5,
                     corvec        = NULL) {

  message_fun("poisthin")

  ## Check Input -------------------------------------------------------------
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nsamp <= nrow(mat))
  assertthat::assert_that(ngene + skip_gene <= ncol(mat))
  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::are_equal(1L,
                        length(skip_gene),
                        length(prop_null),
                        length(alpha),
                        length(group_prop))
  assertthat::assert_that(prop_null >= 0, prop_null <= 1)
  assertthat::assert_that(group_prop >= 0, group_prop <= 1)

  gselect <- match.arg(gselect)

  if (gselect == "custom") {
    stopifnot(is.logical(gvec))
    stopifnot(length(gvec) == ncol(mat))
    stopifnot(sum(gvec) == ngene)
  } else {
    if (!is.null(gvec)) {
      warning('gvec is specified but being ignored since gselect is not "custom"')
    }
  }

  group_assign <- match.arg(group_assign)
  if (group_assign == "cor") {
    stopifnot(!is.null(corvec))
    stopifnot(is.numeric(corvec))
    stopifnot(length(corvec) < min(nsamp, ngene))
    stopifnot(crossprod(corvec) < 1)
  } else {
    stopifnot(is.null(corvec))
  }

  ## get gene indices ---------------------------------------------------------
  med_express <- apply(mat, 2, stats::median)
  mean_express <- colMeans(mat)
  order_vec <- order(med_express, mean_express, decreasing = TRUE)

  if (gselect == "max") {
    gindices <- order_vec[(skip_gene + 1):(skip_gene + ngene)]
  } else if (gselect == "rand_max") {
    first_zero <- match(0, med_express)
    max_gene <- min(c(first_zero, 2 * ngene + skip_gene, ncol(mat)), na.rm = TRUE)
    if (max_gene < ngene + skip_gene) {
      warning("including some low-expressed genes in the sample.")
      max_gene <- ngene + skip_gene
    }
    gindices <- sample(x = order_vec[(skip_gene + 1):max_gene], size = ngene)
  } else if (gselect == "random") {
    gindices <- sample(x = sample(order_vec[(skip_gene + 1):(skip_gene + ngene)]), size = ngene)
  } else if (gselect == "custom") {
    gindices <- (seq_len(ncol(mat)))[gvec]
    if (skip_gene > 0) {
      warning('ignoring skip_gene because gselect = "custom"')
    }
  } else if (gselect == "mean_max") {
    order_vec_means <- order(mean_express, decreasing = TRUE)
    gindices <- order_vec_means[(skip_gene + 1):(skip_gene + ngene)]
  }

  ## Get submat --------------------------------------------------------------
  gindices <- sort(gindices)
  samp_indices <- sort(sample(seq_len(nrow(mat)), size = nsamp))
  submat <- mat[samp_indices, gindices, drop = FALSE]

  ## Group assignment --------------------------------------------------------
  if (group_assign == "frac") {
    group_indicator <- rep(FALSE, length = nsamp)
    group_indicator[sample(seq_len(nsamp), size = round(nsamp * group_prop))] <- TRUE
  } else if (group_assign == "random") {
    group_indicator <- sample(x = c(TRUE, FALSE),
                              size = nsamp,
                              replace = TRUE,
                              prob = c(group_prop, 1 - group_prop))
  } else if (group_assign == "cor") {
    cout <- corassign(mat    = submat,
                      nfac   = length(corvec),
                      corvec = corvec,
                      return = "full")
    group_indicator <- cout$x == 1L
  } else {
    stop("poisthin: how did you get here?")
  }

  if (all(group_indicator) | all(!group_indicator)) {
    warning("All samples were assigned to the same group.")
  }

  ## Draw signal -------------------------------------------------------------
  nsignal <- round(ngene * (1 - prop_null))
  if (nsignal > 0) {
    signal_params$n <- nsignal
    signal_vec      <- do.call(what = signal_fun, args = signal_params) ## log2-fold change

    stopifnot(length(signal_vec) == nsignal)

    which_signal <- sort(sample(seq_len(ncol(submat)), nsignal)) # location of signal

    ## Deal with alpha here ----------------------------
    if (abs(alpha) > 10 ^ -6) {
      sd_vec <- apply(log2(submat[, which_signal, drop = FALSE] + 1), 2, stats::sd) / sqrt(nrow(submat))
      assertthat::are_equal(length(sd_vec), length(signal_vec))
      signal_vec <- signal_vec * (sd_vec ^ alpha)
    }

    sign_vec  <- sign(signal_vec) # sign of signal
    bin_probs <- 2 ^ -abs(signal_vec) # binomial prob

    ng1 <- sum(group_indicator)
    ng2 <- sum(!group_indicator)

    submat[group_indicator, which_signal[sign_vec > 0]] <-
      matrix(stats::rbinom(n    = sum(sign_vec > 0) * ng1,
                           size = c(submat[group_indicator, which_signal[sign_vec > 0]]),
                           prob = rep(bin_probs[sign_vec > 0], each = ng1)),
             nrow = ng1)

    submat[!group_indicator, which_signal[sign_vec < 0]] <-
      matrix(stats::rbinom(n    = sum(sign_vec < 0) * ng2,
                           size = c(submat[!group_indicator, which_signal[sign_vec < 0]]),
                           prob = rep(bin_probs[sign_vec < 0], each = ng2)),
             nrow = ng2)

    beta <- rep(0, ngene)
    beta[which_signal] <- -1 * signal_vec ## -1 because of way design matrix is created
  } else if (nsignal == 0 & abs(prop_null - 1) > 10 ^ -6) {
    warning("no genes were given signal since (1 - prop_null) * ngene was very close to zero")
    beta <- rep(0, ngene)
  } else {
    beta <- rep(0, ngene)
  }

  X <- stats::model.matrix(~group_indicator)
  return_list <- list(Y = submat, X = X, beta = beta)

  if (group_assign == "cor") {
    return_list$corassign <- cout
  }

  return(return_list)
}


#' Group assignment that is correlated with latent factors.
#'
#' We extract latent factors from the log of \code{mat} using an SVD, then
#' generate an underlying group-assignment variable from a conditional
#' normal distribution (conditional on the latent factors). This underlying
#' group-assignment variable is used to assign groups.
#'
#' If \code{nfac} is provided, then \code{corvec} must be the same length as \code{nfac}.
#' If \code{nfac} is not provided, then it is assumed that the first \code{nfac}
#' elements of \code{corvec} are the underlying correlations, if \code{nfac} turns out to be
#' smaller than the length of \code{corvec}. If \code{nfac} turns
#' out to be larger than the length of \code{corvec}, then the factors without
#' defined correlations are assumed to have correlation 0.
#'
#' @param mat A matrix of count data. The rows index the individuals and
#'     the columns index the genes.
#' @param nfac The number of latent factors. If \code{NULL}, then we will
#'     use \code{\link[cate]{est.factor.num}} from the cate
#'     package to choose the number of latent factors.
#' @param corvec The vector of correlations. \code{corvec[i]} is the correlation
#'     between latent factor \code{i} and the underlying group-assignment variable.
#'     You can think of the correlations in \code{corvec} as a kind of "tetrachoric
#'     correlation." If \code{NULL}, then it assumes independence between
#'     factors and group assignment. Note that the correlations of the
#'     latent factors with the observed group-assignment vector (instead of the
#'     latent group-assignment vector) will be \code{corvec * sqrt(2 / pi)}.
#' @param return What should we return? Just the group assignment
#'     (\code{"group"}) or a list of a bunch of things (\code{"full"}).
#'
#' @return A list with some or all of the following elements:
#'     \describe{
#'       \item{\code{x}}{The vector of group assignments. \code{0L} indicates
#'           membership to one group and \code{1L} indicates membership to
#'           the other group.}
#'       \item{\code{nfac}}{The number of assumed latent factors.}
#'       \item{\code{facmat}}{A matrix, whose columns contain the latent factors.}
#'       \item{\code{groupfac}}{The underlying group-assignment factor.}
#'       \item{\code{corvec}}{The correlation vector. Note that this is the
#'           correlation between random variables observed in \code{groupfac}
#'           and \code{facmat}, }
#'     }
#'     If \code{return = "group"}, then the list only contains \code{x}.
#'
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Jingshu Wang and Qingyuan Zhao (2015). cate: High Dimensional
#'     Factor Analysis and Confounder Adjusted Testing and Estimation.
#'     R package version 1.0.4. \url{https://cran.r-project.org/package=cate}}
#' }
#'
#' @examples
#' ## Simulate data from given matrix of counts
#' ## In practice, you would obtain Y from a real dataset, not simulate it.
#' set.seed(1)
#' nsamp <- 1000
#' ngene <- 10
#' Y <- matrix(stats::rpois(nsamp * ngene, lambda = 50), nrow = ngene)
#'
#' ## Set target correlation to be 0.9 and nfac to be 1
#' corvec <- 0.9
#' nfac   <- 1
#'
#' ## Group assignment
#' cout <- corassign(mat    = t(Y),
#'                   nfac   = nfac,
#'                   corvec = corvec,
#'                   return = "full")
#'
#' ## Correlation between facmat and groupfac should be about 0.9
#' cor(cout$facmat, cout$groupfac)
#'
#' ## Correlation between facmat and x should be about 0.9 * sqrt(2 / pi)
#' cor(cout$facmat, cout$x)
#' corvec * sqrt(2 / pi)
#'
#'
#' @export
corassign <- function(mat,
                      nfac   = NULL,
                      corvec = NULL,
                      return = c("group", "full")) {

  ## Check input --------------------------------------------------------------
  return <- match.arg(return)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(all(!is.na(mat)))
  assertthat::assert_that(all(mat >= 0))

  if (!is.null(nfac)) {
    stopifnot(is.numeric(nfac))
    stopifnot(length(nfac) == 1)
    stopifnot(nfac >= 0)
    stopifnot(nfac < min(nrow(mat), ncol(mat)))
  }

  if (!is.null(corvec)) {
    stopifnot(is.numeric(corvec))
    stopifnot(crossprod(corvec) < 1)
  }

  n <- nrow(mat)
  p <- ncol(mat)

  ## Get residuals ------------------------------------------------------------
  resmat <- stats::resid(stats::lm(log2(mat + 1) ~ 1))

  ## Find number of latent factors if not provided ----------------------------
  if (is.null(nfac)) {
    nfac <- cate::est.factor.num(resmat,
                                 method = "bcv",
                                 bcv.plot = FALSE)$r

    ## Pad or delete corvec where necessary.
    if (is.null(corvec)) {
      corvec <- rep_len(x = 0, length.out = nfac)
    } else if (nfac < length(corvec)) {
      corvec <- corvec[seq_len(nfac)]
      message(paste0("Only using first ", nfac, " elements of corvec."))
    } else if (nfac > length(corvec)) {
      npad <- nfac - length(corvec)
      corvec <- c(corvec, rep_len(0, npad))
      message(paste0("Padding last ", npad, " (of ", nfac, " elements) of corvec with 0's."))
    } else {
      message(paste0("Using all ", length(corvec), " elements of corvec."))
    }
  } else {
    if (is.null(corvec)) {
      corvec <- rep_len(x = 0, length.out = nfac)
    } else {
      stopifnot(nfac == length(corvec))
    }
  }

  ## Generate group assignment depending on if nfac == 0 and corvec == 0
  if (nfac == 0) {
    aslist <- uncorassign(n, return = "full")
    x <- aslist$x
    w <- aslist$groupfac
    if (return == "full") {
      facmat <- matrix(nrow = n, ncol = 0)
    }
  } else if (all(corvec == 0)) {
    aslist <- uncorassign(n, return = "full")
    x <- aslist$x
    w <- aslist$groupfac
    if (return == "full") { ## tiny optimization
      facmat <- irlba::irlba(A = resmat, nv = 0, nu = nfac)$u * sqrt(n)
    }
  } else {
    ## Get factors and normalize ----------------------------------------------
    facmat <- irlba::irlba(A = resmat, nv = 0, nu = nfac)$u * sqrt(n)

    ## Generate assignment factor ---------------------------------------------
    w <- stats::rnorm(n    = n,
                      mean = c(facmat %*% corvec),
                      sd   = sqrt(1 - c(crossprod(corvec))))

    ## Group assignment based on assignment factor and return -----------------
    x <- ifelse(w > 0, 1L, 0L)
  }

  ## Return -------------------------------------------------------------------
  return_list   <- list()
  return_list$x <- x
  if (return == "full") {
    return_list$nfac     <- nfac
    return_list$facmat   <- facmat
    return_list$groupfac <- w
    return_list$corvec   <- corvec
  }
  return(return_list)
}


#' Group assignment independent of anything.
#'
#' @param n The sample size.
#' @param return Should we just return a list with just the
#'     vector of assignment (\code{"group"})
#'     or a list with the vector of assignments and the vector of latent
#'     variables (\code{"full"})?
#'
#' @return A list with some or all of the following elements.
#'     \describe{
#'       \item{\code{x}}{The group assignment. \code{1L} for one group and
#'            \code{0L} for the other group.}
#'       \item{\code{w}}{The latent assignment vector (only returned if
#'            \code{return = "full"}). Negative corresponds to one group
#'            and positive corresponds to the other group.}
#'     }
#'
#' @author David Gerard
uncorassign <- function(n,
                        return = c("group", "full")) {
  assertthat::assert_that(n > 0)
  return <- match.arg(return)
  ret_vec   <- list()
  if (return == "group") {
    ret_vec$x <- sample(x       = c(0L, 1L),
                        size    = n,
                        replace = TRUE)
  } else if (return == "full") {
    ret_vec$groupfac <- stats::rnorm(n = n)
    ret_vec$x <- ifelse(ret_vec$groupfac > 0, 1L, 0L)
  }
  return(ret_vec)
}
