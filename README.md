# EBGOAL

**E-Value and Ball Covariance-Based Generalized Outcome-Adaptive LASSO**
  
## Overview
  
  EBGOAL is an R package designed for covariate selection in high-dimensional observational studies, with a particular focus on identifying appropriate covariate adjustment sets when nonlinear relationships may exist between covariates and outcomes. It provides a framework for supporting propensity score analysis and investigating potential causal relationships in complex biomedical data.


## Requirements

EBGOAL depends on the following R packages:
  

install.packages(c(
  "MASS",
  "Matrix",
  "survey",
  "EValue",
  "lqa",
  "Ball",
  "CBPS"
))


## Usage

Prepare a data frame containing the outcome, binary treatment, and candidate covariates. Variables used by the procedure should be complete and numeric. Categorical variables should be converted to numeric indicator variables, and the treatment must be coded as `0` and `1`.


source("R/ATE.R")
source("R/EBGOAL.R")

fit <- EBGOAL(
  Data = analysis_data,
  treatment = "A",
  outcome = "Y",
  var.list = candidate_variables,
  lambda_vec = c(
    -10, -5, -2, -1, -0.75,
    -0.5, -0.25, 0.25, 0.49
  ),
  S_lam = c(
    0,
    10^c(
      -2, -1.5, -1, -0.75, -0.5,
      -0.25, 0, 0.25, 0.5, 1
    )
  ),
  gamma_convergence_factor = -2,
  corprior = 0.01
)

fit$ATE
fit$SE
fit$selectedVar
fit$lambda1
fit$lambda2
fit$bweight


The main outputs include the estimated treatment effect (`ATE`), its standard error (`SE`), selected covariates (`selectedVar`), selected tuning parameters (`lambda1` and `lambda2`), and estimated propensity weights (`bweight`).


## Application

EBGOAL can be applied to high-dimensional biomedical studies, including genomic and other omics research, to explore potential causal relationships and underlying biological mechanisms. Potential applications include identifying mutation-associated genes, investigating their relationships with molecular and clinical phenotypes, and exploring biological pathways involved in disease development and progression.