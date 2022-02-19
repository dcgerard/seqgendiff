####################
## Functions for selecting and filtering a count matrix
####################

#' Subsample the rows and columns of a count matrix.
#'
#' It is a good idea to subsample (each iteration) the genes and samples from
#' a real RNA-seq dataset prior to applying \code{\link{thin_diff}}
#' (and related functions) so that your conclusions are not dependent on the
#' specific structure of your dataset. This function is designed to efficiently
#' do this for you.
#'
#' The samples (columns) are chosen randomly, with each sample having
#' an equal probability of being in the sub-matrix. The genes are selected
#' according to one of four schemes (see the description of the \code{gselect}
#' argument).
#'
#' If you have edgeR installed, then some functionality is provided for
#' filtering out the lowest expressed genes prior to applying subsampling
#' (see the \code{filter_first} argument).
#' This filtering scheme is described in Chen et al. (2016).
#' If you want more control over this filtering, you should use
#' the \code{\link[edgeR]{filterByExpr}} function from edgeR directly. You
#' can install edgeR by following instructions at
#' \doi{10.18129/B9.bioc.edgeR}.
#'
#' @param mat A numeric matrix of RNA-seq counts. The rows index the genes
#'     and the columns index the samples.
#' @param nsamp The number of samples (columns) to select from \code{mat}.
#' @param ngene The number of genes (rows) to select from \code{mat}.
#' @param gselect How should we select the subset of genes? Options include:
#'     \describe{
#'       \item{\code{random}}{Randomly select the genes, with each gene having
#'           an equal probability of being included in the subsampled matrix.}
#'       \item{\code{max}}{Choose the \code{ngene} most median-expressed genes.
#'           Ties are broken by mean-expression.}
#'       \item{\code{mean_max}}{Choose the \code{ngene} most mean-expressed
#'           genes.}
#'       \item{\code{custom}}{A user-specified list of genes. If
#'           \code{gselect = "custom"} then \code{gvec} needs to be
#'           non-\code{NULL}.}
#'     }
#' @param gvec A logical vector of length \code{nrow(mat)}. A \code{TRUE}
#'     in position \eqn{i} indicates inclusion into the smaller dataset.
#'     Hence, \code{sum(gvec)} should equal \code{ngene}.
#' @param filter_first Should we first filter genes by the method of
#'     Chen et al. (2016) (\code{TRUE}) or not (\code{FALSE})? If
#'     \code{TRUE} then the \code{edgeR} package should be installed.
#' @param nskip The number of median-maximally expressed genes to skip.
#'     Not used if \code{gselect = "custom"}.
#'
#' @return A numeric matrix, which is a \code{ngene} by \code{nsamp} sub-matrix
#'     of \code{mat}. If \code{rownames(mat)} is \code{NULL}, then the
#'     row names of the returned matrix are the indices in \code{mat} of the
#'     selected genes. If \code{colnames(mat)} is \code{NULL}, then the
#'     column names of the returned matrix are the indices in \code{mat} of
#'     the selected samples.
#'
#' @references
#' \itemize{
#'   \item{Chen, Yunshun, Aaron TL Lun, and Gordon K. Smyth. "From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline." \emph{F1000Research} 5 (2016). \doi{10.12688/f1000research.8987.2}.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#'
#' @examples
#' ## Simulate data from given matrix of counts
#' ## In practice, you would obtain mat from a real dataset, not simulate it.
#' set.seed(1)
#' n   <- 100
#' p   <- 1000
#' mat <- matrix(stats::rpois(n * p, lambda = 50), nrow = p)
#'
#' ## Subsample the matrix, then feed it into a thinning function
#' submat <- select_counts(mat = mat, nsamp = 10, ngene = 100)
#' thout  <- thin_2group(mat = submat, prop_null = 0.5)
#'
#' ## The rownames and colnames (if NULL in mat) tell you which genes/samples
#' ## were selected.
#' rownames(submat)
#' colnames(submat)
#'
select_counts <- function(mat,
                          nsamp        = ncol(mat),
                          ngene        = nrow(mat),
                          gselect      = c("random",
                                           "max",
                                           "mean_max",
                                           "custom"),
                          gvec         = NULL,
                          filter_first = FALSE,
                          nskip        = 0L) {
  ## Check input --------------------------------------------------------------
  gselect <- match.arg(gselect)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::is.count(nsamp)
  assertthat::is.count(ngene)
  assertthat::is.count(nskip)
  assertthat::assert_that(ncol(mat) >= nsamp)
  assertthat::assert_that(nrow(mat) >= ngene + nskip)
  if (!is.null(gvec)) {
    stopifnot(is.logical(gvec))
    stopifnot(length(gvec) == nrow(mat))
    stopifnot(sum(gvec) == ngene)
    if (gselect != "custom") {
      warning("Ignoring gvec since `gselect != \"custom\"`")
    } else if (nskip != 0L) {
      warning("`gvec` is non-NULL, so ignoring `nskip`")
    }
  }
  assertthat::assert_that(is.logical(filter_first))
  assertthat::are_equal(1L, length(filter_first))
  if (filter_first & gselect == "custom") {
    stop("Cannot have both `filter_first = TRUE` and `gselect = \"custom\"")
  }

  ## Filter mat with edgeR ---------------------------------------------------
  if (filter_first & !requireNamespace("edgeR", quietly = TRUE)) {
    stop(paste0("If `filter_first = TRUE`, then edgeR needs to be installed\n",
                "You can install edgeR with:\n\n",
                "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
                "    install.packages(\"BiocManager\")\n\n",
                "BiocManager::install(\"edgeR\")"))
  } else if (filter_first) {
    eout <- edgeR::filterByExpr(y = mat)
    mat <- mat[eout, , drop = FALSE]
    if (nrow(mat) < ngene + nskip) {
      stop(paste0("Filtering resulting in too few genes\n",
                  "Try decreasing `ngene` or `nskip`\n",
                  "Or try setting `filter_first = FALSE`"))
    }
  }

  ## Select samples -----------------------------------------------------------
  which_samp <- sort(sample(x = seq_len(ncol(mat)), size = nsamp))

  ## Select genes -------------------------------------------------------------
  medvec  <- apply(mat, 1, stats::median, na.rm = TRUE)
  meanvec <- rowMeans(mat, na.rm = TRUE)
  if (nskip == 0L) {
    keep_vec <- seq_len(nrow(mat))
  } else {
    keep_vec <- seq_len(nrow(mat))[-seq_len(nskip)]
  }

  if (gselect == "custom") {
    which_gene <- seq_len(nrow(mat))[gvec]
  } else if (gselect == "max") {
    which_gene <- sort(order(medvec, meanvec, decreasing = TRUE)[keep_vec][seq_len(ngene)])
  } else if (gselect == "mean_max") {
    which_genes <- sort(order(medvec, meanvec, decreasing = TRUE)[keep_vec])
    which_gene  <- sort(which_genes[order(meanvec[which_genes], decreasing = TRUE)][seq_len(ngene)])
  } else if (gselect == "random") {
    which_gene <- order(medvec, meanvec, decreasing = TRUE)[keep_vec]
    which_gene <- sort(sample(which_gene, size = ngene))
  }

  ## Subsample matrix ---------------------------------------------------------
  submat <- mat[which_gene, which_samp, drop = FALSE]

  ## Give rownames and colnames -----------------------------------------------
  if (is.null(rownames(submat))) {
    rownames(submat) <- which_gene
  }
  if (is.null(colnames(submat))) {
    colnames(submat) <- which_samp
  }

  ## Return -------------------------------------------------------------------
  return(submat)
}















