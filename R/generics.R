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
#' This only keeps the \code{mat} and \code{designmat} elements of
#' the ThinData object. Notably, the new SummarizedExperiment object
#' does not contain the true coefficients.
#'
#' @param obj A ThinData S3 object. This is generally output by either
#'     \code{\link{thin_diff}}, \code{\link{thin_2group}}, or
#'     \code{\link{thin_lib}}.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} S4
#'     object. This is often used in Bioconductor when performing
#'     differential expression analysis.
#'
#' @export
#'
#' @author David Gerard
ThinDataToSummarizedExperiment <- function(obj) {
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    assertthat::assert_that(is.ThinData(obj))
    se <- SummarizedExperiment::SummarizedExperiment(assays = obj$mat,
                                                     colData = cbind(1, obj$designmat))
  } else {
    warning(paste0("Need to install SummarizedExperiment to use ThinDataToSummarizedExperiment()\n",
                   "Type in R:\n\n",
                   "install.packages('BiocManager')\n",
                   "BiocManager::install('SummarizedExperiment')\n\n",
                   "Returning NULL for now."))
    se <- NULL
  }
  return(se)
}




