testthat::test_that(
  "EBGOAL selects covariates and estimates treatment effects",
  {
    # Load the fixed test dataset included in the package
    utils::data(
      "EBGOAL_test",
      package = "EBGOAL",
      envir = environment()
    )
    
    dat <- EBGOAL_test
    
    testthat::expect_equal(nrow(dat), 100L)
    testthat::expect_true(
      all(c("A", "Y", paste0("X", 1:30)) %in% names(dat))
    )
    
    fit <- EBGOAL(
      dat,
      treatment = "A",
      outcome = "Y",
      var.list = paste0("X", 1:30)
    )
    
    testthat::expect_type(fit, "list")
    
    testthat::expect_true(
      all(c("selectedVar", "bweight") %in% names(fit))
    )
    
    testthat::expect_length(
      fit$bweight,
      nrow(dat)
    )
    
    testthat::expect_type(
      fit$bweight,
      "double"
    )
    
    testthat::expect_true(
      all(is.finite(fit$bweight))
    )
    
    dat$bweight <- fit$bweight
    
    ate <- ATE_est(
      dat,
      w = "bweight",
      out = "Y",
      trt.list = "A"
    )
    
    testthat::expect_false(
      is.null(ate$Trt_coef)
    )
  }
)
