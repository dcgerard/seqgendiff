## Generate counts from the theoretical model.

#' Simulate RNA seq data via a factor-augmented regression model on the log-counts.
#'
#' This function simulates RNA-seq data generated under the assumed
#' model on the log-counts (\eqn{Y}):
#' \deqn{Y = XB + WC + ZA + E,}
#' where \eqn{Y} is a matrix of gene-expression levels
#' (rows index samples and columns index genes),
#' \eqn{X} is a matrix of observed covariates of interest,
#' \eqn{B} is a matrix of unobserved coefficients of interest,
#' \eqn{W} is a matrix of observed nuisance covariates,
#' \eqn{C} is a matrix of unobserved coefficients of \eqn{W},
#' \eqn{Z} is a matrix of unobserved confounders, \eqn{A} is a matrix of unobserved
#' coefficients of \eqn{Z}, and \eqn{E} is Gaussian noise with column-specific variances.
#' Various options are available for correlation between X and Z and the proportion of variance
#' explained by X and Z. X and Z can consist of either normal variates or columns of indicator
#' variables. Abilities to vary library size and the effects of transcript length (special cases
#' of unobserved confounding) are included.
#'
#' @author David Gerard
#'
#' @export
#'
#' @param nsamp The sample size.
#' @param ngene The number of genes.
#' @param signal_fun A function that returns a vector of the signal of interest. It must
#'     take as a parameter \code{n}, for the number of samples to return.
#' @param signal_params A list of parameters to pass to \code{signal_fun}.
#'     It cannot include a parameter named \code{n}.
#' @param prop_null The proportion of genes that are null.
#'
theogen <- function(nsamp, ngene, signal_fun = stats::rnorm,
                    signal_params = list(mean = 0, sd = 1),
                    prop_null = 1 ) {



}
