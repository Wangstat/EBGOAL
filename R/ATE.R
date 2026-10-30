#' Estimate a Weighted Treatment Effect Using IPTW
#' @param data Data frame containing the outcome, treatment variables, and weights.
#' @param w Character string naming the weight column.
#' @param out Character string naming the numeric outcome variable.
#' @param trt.list Character vector naming the treatment variables.
#' @return A list containing \code{Trt_coef}, the coefficient matrix from the
#'   survey-weighted model summary.
#' @export
ATE_est <- function(data, w, out, trt.list) {
  
  data <- as.data.frame(data)
  
  design.b <- survey::svydesign(
    ids = ~1,
    weights = stats::reformulate(w),
    data = data
  )
  
  ATE.formula <- stats::reformulate(
    trt.list,
    response = out
  )
  
  tt_ATE <- survey::svyglm(
    ATE.formula,
    design = design.b
  )
  
  list(
    Trt_coef = summary(tt_ATE)$coefficients
  )
}