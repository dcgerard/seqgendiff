#####################
## Methods for the ThinData S3 class.
#####################


is.ThinData <- function(x) {
  inherits(x, "ThinData")
}


#' Provide summary output of a ThinData S3 object.
#'
#' @param object A ThinData S3 object. This is generally output by either
#'     \code{\link{thin_diff}}, \code{\link{thin_2group}}, or
#'     \code{\link{thin_lib}}.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @return Returns nothing. Prints out some summary information on
#'     \code{object}.
#'
#' @export
summary.ThinData <- function(object, ...) {
cat(
  "class: ThinData\n",
  "mat\n",
  "   modified count matrix\n",
  "   dim:", dim(object$mat), "\n",
  "designmat\n",
  "   design matrix\n",
  "   dim:", dim(object$designmat), "\n",
  "coefmat\n",
  "   coefficients of design matrix\n",
  "   dim:", dim(object$coefmat), "\n",
  "design_obs\n",
  "   additional variables to include in the design matrix\n",
  "   dim:", dim(object$design_obs), "\n",
  "sv\n",
  "   estimated surrogate variable(s)\n",
  "   dim:", dim(object$sv), "\n",
  "cormat\n",
  "   theoretical correlation(s) between sv and matching_var\n",
  "   rows = observed covariates, columns = surrogate variables\n",
  "   dim:", dim(object$cormat), "\n",
  "matching_var\n",
  "   simulated variable(s) used to permute elements of designmat\n",
  "   dim:", dim(object$matching_var), "\n"
)
}


#' Converts a ThinData S3 object into a SummarizedExperiment S4 object.
#'
#' This only keeps the \code{mat}, \code{design_obs}, \code{designmat},
#' and \code{coefmat} elements of the ThinData object.
#'
#' @param obj A ThinData S3 object. This is generally output by either
#'     \code{\link{thin_diff}}, \code{\link{thin_2group}},
#'     \code{\link{thin_lib}}, \code{\link{thin_gene}}, or
#'     \code{\link{thin_all}}.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} S4
#'     object. This is often used in Bioconductor when performing
#'     differential expression analysis.
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#' ## Generate simulated data and modify using thin_diff().
#' ## In practice, you would use real data, not simulated.
#' set.seed(1)
#' n <- 10
#' p <- 1000
#' Z <- matrix(abs(rnorm(n, sd = 4)))
#' alpha <- matrix(abs(rnorm(p, sd = 1)))
#' mat <- round(2^(alpha %*% t(Z) + abs(matrix(rnorm(n * p, sd = 5),
#'                                             nrow = p,
#'                                             ncol = n))))
#' design_perm <- cbind(rep(c(0, 1), length.out = n), runif(n))
#' coef_perm   <- matrix(rnorm(p * ncol(design_perm), sd = 6), nrow = p)
#' design_obs  <- matrix(rnorm(n), ncol = 1)
#' target_cor <- matrix(c(0.9, 0))
#' thout <- thin_diff(mat            = mat,
#'                    design_perm    = design_perm,
#'                    coef_perm      = coef_perm,
#'                    target_cor     = target_cor,
#'                    design_obs     = design_obs,
#'                    permute_method = "hungarian")
#'
#' ## Convert ThinData object to SummarizedExperiment object.
#' seobj <- ThinDataToSummarizedExperiment(thout)
#' class(seobj)
#'
#' ## The "O1" variable in the colData corresponds to design_obs.
#' ## The "P1" and "P2" variables in colData correspond to design_perm.
#' seobj
#' }
#'
#' @author David Gerard
ThinDataToSummarizedExperiment <- function(obj) {
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    assertthat::assert_that(is.ThinData(obj))
    if (ncol(obj$coefmat) > 0) {
      colnames(obj$coefmat) <- paste0("true_", colnames(obj$coefmat))
    }
    rownames(obj$coefmat) <- paste0("gene", seq_len(nrow(obj$coefmat)))
    rownames(obj$mat)     <- paste0("gene", seq_len(nrow(obj$mat)))
    ## drop intercept
    overall_design <- cbind(obj$design_obs[, -1, drop = FALSE], obj$designmat)
    rownames(overall_design) <- paste0("sample", seq_len(nrow(overall_design)))
    colnames(obj$mat) <- paste0("sample", seq_len(ncol(obj$mat)))
    se <- SummarizedExperiment::SummarizedExperiment(assays = obj$mat,
                                                     colData = overall_design,
                                                     rowData = obj$coefmat)
  } else {
    warning(paste0("Need to install SummarizedExperiment to use ",
                   "ThinDataToSummarizedExperiment()\n",
                   "Type in R:\n\n",
                   "install.packages('BiocManager')\n",
                   "BiocManager::install('SummarizedExperiment')\n\n",
                   "Returning NULL for now."))
    se <- NULL
  }
  return(se)
}

#' Converts a ThinData S3 object into a DESeqDataSet S4 object.
#'
#' The design formula in the resulting DESeqDataSet is just the sum of all
#' variables in \code{designmat} from the ThinData object (except the
#' intercept term). You should change this design formula if you want to
#' study other models.
#'
#' @inheritParams ThinDataToSummarizedExperiment
#'
#' @return A \code{\link[DESeq2]{DESeqDataSet}} S4
#'     object. This will allow you to insert the simulated
#'     data directly into DESeq2.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Generate simulated data and modify using thin_diff().
#' ## In practice, you would use real data, not simulated.
#' set.seed(1)
#' n <- 10
#' p <- 1000
#' Z <- matrix(abs(rnorm(n, sd = 4)))
#' alpha <- matrix(abs(rnorm(p, sd = 1)))
#' mat <- round(2^(alpha %*% t(Z) + abs(matrix(rnorm(n * p, sd = 5),
#'                                             nrow = p,
#'                                             ncol = n))))
#' design_perm <- cbind(rep(c(0, 1), length.out = n), runif(n))
#' coef_perm   <- matrix(rnorm(p * ncol(design_perm), sd = 6), nrow = p)
#' design_obs  <- matrix(rnorm(n), ncol = 1)
#' target_cor <- matrix(c(0.9, 0))
#' thout <- thin_diff(mat            = mat,
#'                    design_perm    = design_perm,
#'                    coef_perm      = coef_perm,
#'                    target_cor     = target_cor,
#'                    design_obs     = design_obs,
#'                    permute_method = "hungarian")
#'
#' ## Convert ThinData object to DESeqDataSet object.
#' seobj <- ThinDataToDESeqDataSet(thout)
#' class(seobj)
#'
#' ## The "O1" variable in the colData corresponds to design_obs.
#' ## The "P1" and "P2" variables in colData correspond to design_perm.
#' seobj
#' }
#'
#' @author David Gerard
#'
ThinDataToDESeqDataSet <- function(obj) {
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    se <- ThinDataToSummarizedExperiment(obj)
    design_form <- stats::formula(paste0("~ ", paste0(paste0("`", names(SummarizedExperiment::colData(se)), "`"), collapse = " + ")))
    dds <- DESeq2::DESeqDataSet(se = se, design = design_form)
  } else {
    warning(paste0("Need to install DESeq2 to use ThinDataToDESeqDataSet()\n",
                   "Type in R:\n\n",
                   "install.packages('BiocManager')\n",
                   "BiocManager::install('DESeq2')\n\n",
                   "Returning NULL for now."))
    dds <- NULL
  }
  return(dds)
}
