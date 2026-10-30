#' Simulated dataset for EBGOAL
#'
#' A simulated dataset used to demonstrate and test
#' the EBGOAL high-dimensional propensity-score procedure.
#'
#' @format A data frame with 100 observations and 32 variables:
#' \describe{
#'   \item{A}{Binary treatment indicator (0 or 1).}
#'   \item{Y}{Continuous outcome variable.}
#'   \item{X1-X30}{Thirty covariates used in the EBGOAL procedure.}
#' }
#'
#' @details
#' This artificial dataset is provided for demonstrating
#' the use of EBGOAL.
#'
#' @examples
#' data("EBGOAL_test")
#' dim(EBGOAL_test)
#'
#' @source Simulated data generated for the EBGOAL package.
#'
#' @name EBGOAL_test
#' @docType data
NULL