#' Fit the EBGOAL high-dimensional propensity-score procedure
#'
#' Select covariates using leave-one-out E-value changes and conditional
#' Ball covariance, tune the adaptive-LASSO penalty using weighted mean
#' absolute differences, and estimate propensity weights using npCBPS
#' for treatment effect estimation.
#'
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats sd
#' @importFrom lqa lqa.formula
#' @param Data A complete-case data frame containing the outcome, binary
#'   treatment, and candidate covariates.
#' @param treatment Name of the binary treatment variable coded 0/1.
#' @param outcome Name of the numeric outcome variable.
#' @param var.list Character vector of candidate covariate names.
#' @param lambda_vec Candidate values for the first penalty parameter.
#' @param S_lam Candidate values for the second penalty parameter.
#' @param gamma_convergence_factor Numeric constant used to calculate gamma.
#' @param corprior Correlation-prior parameter passed to npCBPS.
#' @return A list containing ATE, SE, lambda1, lambda2, selectedVar,
#'   and bweight.
#' @export
EBGOAL <- function(
    Data,
    treatment,
    outcome,
    var.list,
    lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
    S_lam = c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1)),
    gamma_convergence_factor = -2,
    corprior = 0.01
) {
  
  start_time <- Sys.time()
  
  calculate_Evalue <- function(ATE, SE, Y) {
    if (!is.finite(ATE) || !is.finite(SE) || SE <= 0 ||
        !is.finite(stats::sd(Y)) || stats::sd(Y) <= 0) {
      return(NA_real_)
    }
    out <- tryCatch(
      EValue::evalues.OLS(
        est = ATE,
        se = SE,
        sd = stats::sd(Y),
        delta = 1
      ),
      error = function(e) NULL
    )
    if (is.null(out) ||
        !("E-values" %in% rownames(out)) ||
        !("point" %in% colnames(out))) {
      return(NA_real_)
    }
    as.numeric(out["E-values", "point"])
  }
  
  Causal.cor <- function(x, y, z, distance = FALSE) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.numeric(z)
    i0 <- which(z == 0)
    i1 <- which(z == 1)
    if (length(i0) < 2L || length(i1) < 2L) {
      return(NA_real_)
    }
    alpha <- length(i0) / length(z)
    if (distance) {
      x0 <- x[i0, i0, drop = FALSE]
      y0 <- y[i0, i0, drop = FALSE]
      x1 <- x[i1, i1, drop = FALSE]
      y1 <- y[i1, i1, drop = FALSE]
      val <- alpha * Ball::bcov(
        x0, y0, distance = TRUE
      )^2 +
        (1 - alpha) * Ball::bcov(
          x1, y1, distance = TRUE
        )^2
    } else {
      x0 <- x[i0, , drop = FALSE]
      y0 <- y[i0, , drop = FALSE]
      x1 <- x[i1, , drop = FALSE]
      y1 <- y[i1, , drop = FALSE]
      val <- alpha * Ball::bcov(x0, y0)^2 +
        (1 - alpha) * Ball::bcov(x1, y1)^2
    }
    as.numeric(val)
  }
  
  wAMD_function <- function(DataM, varlist, trt.var, wgt, beta) {
    trt <- untrt <- diff_vec <- rep(
      NA_real_, length(varlist)
    )
    names(trt) <- names(untrt) <- names(diff_vec) <- varlist
    for (j in seq_along(varlist)) {
      x <- DataM[[varlist[j]]]
      w <- DataM[[wgt]]
      A <- DataM[[trt.var]]
      d1 <- sum(w[A == 1], na.rm = TRUE)
      d0 <- sum(w[A == 0], na.rm = TRUE)
      if (!is.finite(d1) || !is.finite(d0) ||
          d1 <= 0 || d0 <= 0) {
        return(list(
          diff_vec = diff_vec,
          wdiff_vec = rep(NA_real_, length(varlist)),
          wAMD = NA_real_
        ))
      }
      trt[j] <- sum(
        x[A == 1] * w[A == 1],
        na.rm = TRUE
      ) / d1
      untrt[j] <- sum(
        x[A == 0] * w[A == 0],
        na.rm = TRUE
      ) / d0
      diff_vec[j] <- abs(trt[j] - untrt[j])
    }
    wdiff_vec <- diff_vec * abs(beta)
    score <- sum(wdiff_vec)
    list(
      diff_vec = diff_vec,
      wdiff_vec = wdiff_vec,
      wAMD = if (is.finite(score)) score else NA_real_
    )
  }
  
  standardize_matrix <- function(Data, vars) {
    means <- vapply(
      Data[, vars, drop = FALSE],
      mean,
      numeric(1)
    )
    sds <- vapply(
      Data[, vars, drop = FALSE],
      stats::sd,
      numeric(1)
    )
    bad <- which(!is.finite(sds) | sds == 0)
    if (length(bad)) {
      stop(
        "Candidate covariates with zero or invalid SD: ",
        paste(vars[bad], collapse = ", "),
        call. = FALSE
      )
    }
    for (j in seq_along(vars)) {
      Data[[vars[j]]] <-
        (Data[[vars[j]]] - means[j]) / sds[j]
    }
    list(
      Data = Data,
      means = means,
      sds = sds
    )
  }
  
  make_safe_formula <- function(response, predictors) {
    stats::reformulate(
      predictors,
      response = response
    )
  }
  
  Data <- as.data.frame(Data)
  n <- nrow(Data)
  p <- length(var.list)
  
  Data[[treatment]] <- as.numeric(Data[[treatment]])
  Data[[outcome]] <- as.numeric(Data[[outcome]])
  
  cat("EBGOAL started.\n")
  cat("Sample size:", n, "| Candidate covariates:", p, "\n")
  
  cat("Standardizing data...\n")
  standard_result <- standardize_matrix(Data, var.list)
  Data <- standard_result$Data
  
  A <- Data[[treatment]]
  Y <- Data[[outcome]]
  
  cat("Calculating full-model E-value...\n")
  lm.Y <- lm(
    make_safe_formula(outcome, c(treatment, var.list)),
    data = Data
  )
  
  beta_A_full <- coef(lm.Y)[treatment]
  se_A_full <- summary(lm.Y)$coefficients[
    treatment, "Std. Error"
  ]
  
  E_value_ini <- tryCatch({
    tmp <- EValue::evalues.OLS(
      est = beta_A_full,
      se = se_A_full,
      sd = sd(Y),
      delta = 1
    )
    as.numeric(tmp["E-values", "point"])
  }, error = function(e) NA_real_)
  
  cat("Calculating leave-one-out E-values...\n")
  E_value_x <- rep(NA_real_, p)
  names(E_value_x) <- var.list
  
  for (k in seq_along(var.list)) {
    
    fit_k <- tryCatch(
      lm(
        make_safe_formula(
          outcome,
          c(treatment, var.list[-k])
        ),
        data = Data
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_k) ||
        !(treatment %in% names(coef(fit_k)))) {
      next
    }
    
    se_A_k <- tryCatch(
      summary(fit_k)$coefficients[treatment, "Std. Error"],
      error = function(e) NA_real_
    )
    
    E_value_x[k] <- tryCatch({
      if (!is.finite(se_A_k) || se_A_k <= 0) {
        NA_real_
      } else {
        tmp <- EValue::evalues.OLS(
          est = unname(coef(fit_k)[treatment]),
          se = se_A_k,
          sd = sd(Y),
          delta = 1
        )
        as.numeric(tmp["E-values", "point"])
      }
    }, error = function(e) NA_real_)
  }
  
  delta_abs <- abs(E_value_x - E_value_ini)
  max_delta <- suppressWarnings(max(delta_abs, na.rm = TRUE))
  
  delta_norm <- if (!is.finite(max_delta) || max_delta == 0) {
    rep(0, p)
  } else {
    delta_abs / max_delta
  }
  
  names(delta_norm) <- var.list
  
  cat("Calculating conditional Ball covariance...\n")
  ballcor <- rep(NA_real_, p)
  names(ballcor) <- var.list
  
  for (k in seq_along(var.list)) {
    ballcor[k] <- tryCatch(
      Causal.cor(Data[[var.list[k]]], Y, A),
      error = function(e) NA_real_
    )
  }
  
  ball_abs <- abs(ballcor)
  max_ball <- suppressWarnings(max(ball_abs, na.rm = TRUE))
  
  ball_norm <- if (!is.finite(max_ball) || max_ball == 0) {
    rep(0, p)
  } else {
    ball_abs / max_ball
  }
  
  names(ball_norm) <- var.list
  weight_new <- delta_norm * ball_norm
  names(weight_new) <- var.list
  
  gamma_vals <- 2 * (
    gamma_convergence_factor - lambda_vec + 1
  )
  names(gamma_vals) <- as.character(lambda_vec)
  
  wAMD_matrix <- matrix(
    NA_real_,
    nrow = length(S_lam),
    ncol = length(lambda_vec),
    dimnames = list(
      paste0("lambda2_", S_lam),
      paste0("lambda1_", lambda_vec)
    )
  )
  
  selected_list <- list()
  X <- Data[, var.list, drop = FALSE]
  
  cat("Tuning penalty parameters...\n")
  tuning_start <- Sys.time()
  
  for (m in seq_along(S_lam)) {
    
    lambda2 <- S_lam[m]
    Ip <- sqrt(lambda2) * diag(p)
    
    Anp <- c(A, rep(0, p))
    Xnp <- matrix(0, n + p, p)
    
    for (j in seq_len(p)) {
      Xnp[, j] <- c(X[, j], Ip[, j])
    }
    
    newData <- as.data.frame(Xnp)
    names(newData) <- var.list
    newData[[treatment]] <- Anp
    
    w.form <- make_safe_formula(treatment, var.list)
    
    for (lil in seq_along(lambda_vec)) {
      
      lambda1 <- lambda_vec[lil]
      gamma_now <- gamma_vals[as.character(lambda1)]
      
      key <- paste0(
        "lambda2_", lambda2,
        "_lambda1_", lambda1
      )
      
      if (!is.finite(gamma_now) || gamma_now <= 1) {
        next
      }
      
      penalty_weight <- pmax(abs(weight_new), 1e-8)^(-gamma_now)
      
      oal_pen <- lqa::adaptive.lasso(
        lambda = (n + p)^lambda1,
        al.weights = penalty_weight
      )
      
      logit_oal <- tryCatch(
        lqa.formula(
          w.form,
          data = newData,
          penalty = oal_pen,
          family = binomial(logit)
        ),
        error = function(e) NULL
      )
      
      if (is.null(logit_oal)) {
        next
      }
      
      coef_now <- coef(logit_oal)
      coef_gene <- coef_now[var.list]
      coef_gene[!is.finite(coef_gene)] <- 0
      
      selected_now <- var.list[
        abs((1 + lambda2) * coef_gene) > 1e-12
      ]
      
      selected_list[[key]] <- selected_now
      
      Data_cbps <- Data
      Data_cbps[[treatment]] <- factor(
        Data_cbps[[treatment]],
        levels = c(0, 1)
      )
      
      w.model <- make_safe_formula(treatment, selected_now)
      
      fit_cbps <- tryCatch(
        CBPS::npCBPS(
          w.model,
          data = Data_cbps,
          corprior = corprior,
          print.level = 0
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_cbps)) {
        next
      }
      
      weight_now <- as.numeric(fit_cbps$weights)
      
      if (length(weight_now) != n ||
          any(!is.finite(weight_now)) ||
          any(weight_now <= 0)) {
        next
      }
      
      Data$.__TEMP_WEIGHT__ <- weight_now
      
      wAMD_now <- tryCatch(
        wAMD_function(
          DataM = Data,
          varlist = var.list,
          trt.var = treatment,
          wgt = ".__TEMP_WEIGHT__",
          beta = penalty_weight
        ),
        error = function(e) NULL
      )
      
      Data$.__TEMP_WEIGHT__ <- NULL
      
      if (is.null(wAMD_now) || !is.finite(wAMD_now$wAMD)) {
        next
      }
      
      wAMD_matrix[m, lil] <- wAMD_now$wAMD
    }
    
    cat(
      sprintf(
        "Tuning progress: %d/%d lambda2 values completed.\n",
        m, length(S_lam)
      )
    )
  }
  
  cat(
    sprintf(
      "Penalty tuning completed in %.2f minutes.\n",
      as.numeric(difftime(Sys.time(), tuning_start, units = "mins"))
    )
  )
  
  valid <- which(is.finite(wAMD_matrix), arr.ind = TRUE)
  
  if (nrow(valid) == 0) {
    stop("No valid penalty combination was obtained.")
  }
  
  best_index <- valid[
    which.min(wAMD_matrix[valid]),
    ,
    drop = FALSE
  ]
  
  best_m <- best_index[1, 1]
  best_l <- best_index[1, 2]
  
  lambda2_opt <- S_lam[best_m]
  lambda1_opt <- lambda_vec[best_l]
  gamma_opt <- gamma_vals[as.character(lambda1_opt)]
  
  cat(
    "Optimal parameters:",
    "lambda1 =", lambda1_opt,
    "| lambda2 =", lambda2_opt, "\n"
  )
  
  Ip <- sqrt(lambda2_opt) * diag(p)
  Anp <- c(A, rep(0, p))
  Xnp <- matrix(0, n + p, p)
  
  for (j in seq_len(p)) {
    Xnp[, j] <- c(X[, j], Ip[, j])
  }
  
  newData <- as.data.frame(Xnp)
  names(newData) <- var.list
  newData[[treatment]] <- Anp
  
  penalty_weight_opt <- pmax(abs(weight_new), 1e-8)^(-gamma_opt)
  
  final_penalty <- lqa::adaptive.lasso(
    lambda = (n + p)^lambda1_opt,
    al.weights = penalty_weight_opt
  )
  
  final_lqa <- lqa.formula(
    make_safe_formula(treatment, var.list),
    data = newData,
    penalty = final_penalty,
    family = binomial(logit)
  )
  
  final_coef <- coef(final_lqa)
  final_covariate_coef <- final_coef[var.list]
  
  cat("\n===== EBGOAL coefficient diagnostics =====\n")
  
  cat("Number of candidate covariates:", length(var.list), "\n")
  cat("Number of coefficients returned:", length(final_coef), "\n")
  cat("Number of matched coefficients:", sum(!is.na(final_covariate_coef)), "\n")
  cat("Number of missing coefficients:", sum(is.na(final_covariate_coef)), "\n")
  
  cat("\nCoefficient names:\n")
  print(names(final_coef))
  
  cat("\nExtracted gene coefficients:\n")
  print(final_covariate_coef)
  
  cat("\nCoefficient summary:\n")
  print(summary(as.numeric(final_covariate_coef)))
  
  cat("\nNumber of exact zero coefficients:",
      sum(final_covariate_coef == 0, na.rm = TRUE), "\n")
  
  cat("Number of coefficients with abs(value) <= 1e-12:",
      sum(abs(final_covariate_coef) <= 1e-12, na.rm = TRUE), "\n")
  
  cat("lambda2_opt:", lambda2_opt, "\n")
  
  final_covariate_coef[!is.finite(final_covariate_coef)] <- 0
  
  selectedVar <- var.list[
    abs((1 + lambda2_opt) * final_covariate_coef) > 1e-12
  ]
  
  cat("Selected covariates after filtering:", length(selectedVar), "\n")
  cat("Selected variable names:\n")
  print(selectedVar)
  
  cat("==========================================\n")
  
  cat("Selected covariates:", length(selectedVar), "\n")
  cat("Estimating final propensity weights...\n")
  
  Data_cbps <- Data
  Data_cbps[[treatment]] <- factor(
    Data_cbps[[treatment]],
    levels = c(0, 1)
  )
  
  final_cbps <- CBPS::npCBPS(
    make_safe_formula(treatment, selectedVar),
    data = Data_cbps,
    corprior = corprior,
    print.level = 0
  )
  
  bweight <- as.numeric(final_cbps$weights)
  
  cat("Estimating treatment effect...\n")
  
  Data$.ebgoal_weight <- bweight
  
  ATE_result <- ATE_est(
    data = Data,
    w = ".ebgoal_weight",
    out = outcome,
    trt.list = treatment
  )
  
  ATE <- unname(ATE_result$Trt_coef[treatment, "Estimate"])
  SE <- unname(ATE_result$Trt_coef[treatment, "Std. Error"])
  
  cat(
    sprintf(
      "EBGOAL completed in %.2f minutes.\n",
      as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    )
  )
  
  return(list(
    ATE = ATE,
    SE = SE,
    lambda1 = lambda1_opt,
    lambda2 = lambda2_opt,
    selectedVar = selectedVar,
    bweight = bweight
  ))
}