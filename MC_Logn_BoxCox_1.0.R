# ================================================
# Full Monte Carlo Simulation for Box-Cox Cure Model
# ================================================
rm(list=ls())
options(encoding = "UTF-8")
library(copula)
library(dplyr)
library(tidyverse)
library(numDeriv)

logit <- function(p) log(p / (1 - p))
generate_logn_boxcox_data <- function(n = 2000, seed = 2000, rho = 0.12, alpha_true = 1.0, cens_rate = 0.001) {
  set.seed(seed)
  
  # Generate correlated covariates through Gaussian
  normal_cop <- normalCopula(param = rho, dim = 4, dispstr = "ex")
  U <- rCopula(n, normal_cop)
  U <- as.data.frame(U)
  colnames(U) <- paste0("U", 1:4)
  
  X <- U |>
    dplyr::mutate(
      X1 = case_when(
        U1 <= 0.04 ~ "A",
        U1 <= 0.24 ~ "B",
        U1 <= 0.86 ~ "C",
        TRUE ~ "D"
      ),
      X2 = pmin(pmax(U2 * 0.36, 0.04), 0.36),
      X3 = pmax(pmin(qlnorm(U3, meanlog = -1.555, sdlog = 0.58), 0.99), 0.01),
      X4 = pmin(pmax(rpois(n, lambda = 15*U4), 0),22),
      X5 = pmin(pmax(U5 * 2.5, 0), 2.5)
    ) |>
    dplyr::mutate(
      X1 = factor(X1, levels = c("A", "B", "C", "D"))
    )
  
  X_mat <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = X)
  
  # === Tuned Parameters for ~45% cure rate, ~50% censoring rate ===
  beta_true <- c(
    `(Intercept)` = 1.6,    
    X1B = -0.7,
    X1C = -0.49,
    X1D = -0.41,
    X2 = -5.65,        
    X3 = -0.55,
    X4 = -0.06,
    X5 = 0.4
  )
  
  meanlog_true <- 3.2
  sdlog_true <- 1.3
  
  # Compute phi and cure probability
  eta <- as.vector(X_mat %*% beta_true)
  
  if (alpha_true < 1e-5) {
    phi <- exp(eta)
    cure_prob <- exp(-phi)
  } else {
    phi <- exp(eta) / (1 + alpha_true * exp(eta))
    cure_prob <- (1 - alpha_true * phi)^(1 / alpha_true)
  }
  
  cure_prob <- pmax(pmin(cure_prob, 1), 1e-8)
  
  # Generate censoring times from exponential
  C <- rexp(n, rate = cens_rate)  
  
  time_obs <- numeric(n)
  status <- numeric(n)
  cured_true <- numeric(n)
  
  for (i in 1:n) {
    p0i <- cure_prob[i]
    U_star <- runif(1)
    
    if (U_star <= p0i) {
      time_obs[i] <- C[i]
      status[i] <- 0
      cured_true[i] <- 1
    } else {
      U_double_star <- runif(1)
      
      if (alpha_true < 1e-5) {
        arg <- -log(p0i + (1 - p0i) * U_double_star) / phi[i]
      } else {
        arg <- (1 - (p0i + (1 - p0i) * U_double_star)^alpha_true) / (alpha_true * phi[i])
      }
      
      log_t <- meanlog_true + sdlog_true * qnorm(arg)
      t_event <- exp(log_t)
      
      time_obs[i] <- min(t_event, C[i])
      status[i] <- as.numeric(t_event <= C[i])
      cured_true[i] <- 0
    }
  }
  
  df <- data.frame(
    time = time_obs,
    status = status,
    X1 = X$X1,
    X2 = X$X2,
    X3 = X$X3,
    X4 = X$X4,
    X5 = X$X5,
    cured_true = cured_true
  )

  return(df)
}

# ================== Box-Cox NLL ==================
boxcox_nll <- function(par, time, status, X) {
  alpha_raw <- par[1]
  alpha <- plogis(alpha_raw) 
  
  beta <- par[2:(ncol(X)+1)]
  eta <- as.vector(X %*% beta)
  
  meanlog <- par[ncol(X)+2]
  sdlog <- exp(par[length(par)])
  
  F_t <- plnorm(time, meanlog = meanlog, sdlog = sdlog)
  f_t <- dlnorm(time, meanlog = meanlog, sdlog = sdlog)
  
  if (alpha < 1e-5) {
    phi <- exp(eta)
    S_p <- exp(-phi * F_t)
    f_p <- phi * f_t * exp(-phi * F_t)
  } else {
    phi <- exp(eta) / (1 + alpha * exp(eta))
    term <- 1 - alpha * phi * F_t
    term <- pmax(term, 1e-12)
    S_p <- term^(1/alpha)
    f_p <- phi * f_t * term^(1/alpha - 1)
  }
  
  S_p <- pmax(S_p, 1e-12)
  f_p <- pmax(f_p, 1e-12)
  
  ll <- sum(status * log(f_p) + (1 - status) * log(S_p))
  return(-ll)
}

# ================================================
# Monte Carlo Simulation for Box-Cox
# ================================================

monte_carlo_boxcox <- function(n_sim = 1000, n = 2000, rho = 0.12, alpha_true=1.0) {
  
  results_list <- list()
  
  for (i in 1:n_sim) {
    cat("Simulation", i, "of", n_sim, "\r")
    
    # === Generate data from Box-Cox ===
    dat <- generate_logn_boxcox_data(n = n, seed = 1000 + i,alpha_true=alpha_true, rho = rho)
    
    X <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = dat)
    
    start <- c(qlogis(0.5), rep(0, ncol(X)), 3.5, 1.25)
    
    fit <- optim(par = start, fn = boxcox_nll,
                 time = dat$time, status = dat$status, X = X,
                 method = "BFGS", 
                 control = list(maxit = 10000, reltol = 1e-10),
                 hessian = TRUE)
    
    est <- fit$par
    alpha_hat <- plogis(est[1])
    
     H <- hessian(boxcox_nll,
             fit$par,
             time=dat$time,
             status=dat$status,
             X=X)
     V <- tryCatch(solve(H), error = function(e) NULL)
    
   # V <- tryCatch(solve(fit$hessian), error = function(e) NULL)
    se <- if(!is.null(V)) sqrt(diag(V)) else rep(NA, length(est))
    se_alpha <- alpha_hat * (1 - alpha_hat) * se[1]
    
    results_list[[i]] <- c(alpha = alpha_hat, est[-1], se_alpha = se_alpha, se[-1])
  }
  
  # Parameter names
  param_names <- c("alpha", colnames(X), "meanlog", "log_sdlog")
  
  res_mat <- do.call(rbind, results_list)
  colnames(res_mat) <- c(paste0(param_names, "_est"), paste0(param_names, "_se"))
  res_df <- as.data.frame(res_mat)
  
  # True values from data generation
  true_values <- c(
    alpha = 1.0,
    `(Intercept)` = 1.6,    
    X1B = -0.7,
    X1C = -0.49,
    X1D = -0.41,
    X2 = -5.65,        
    X3 = -0.55,
    X4 = -0.06,
    X5 = 0.4,
    meanlog = 3.2,
    log_sdlog = log(1.3)
  )
  
  # Performance metrics
  metrics <- data.frame(
    Parameter = param_names,
    True      = true_values,
    MeanEst   = colMeans(res_df[, paste0(param_names, "_est")], na.rm = TRUE),
    Bias      = colMeans(res_df[, paste0(param_names, "_est")], na.rm = TRUE) - true_values,
    PctBias   = 100 * (colMeans(res_df[, paste0(param_names, "_est")], na.rm = TRUE) - true_values) / abs(true_values),
    MSE       = colMeans((res_df[, paste0(param_names, "_est")] - true_values)^2, na.rm = TRUE),
    RMSE      = sqrt(colMeans((res_df[, paste0(param_names, "_est")] - true_values)^2, na.rm = TRUE))
  )
  
  # Coverage probabilities
  metrics$CP90 <- NA_real_
  metrics$CP95 <- NA_real_
  
  for (p in param_names) {
    est_col <- paste0(p, "_est")
    se_col  <- paste0(p, "_se")
    metrics$CP90[metrics$Parameter == p] <- mean(
      res_df[[est_col]] - 1.645 * res_df[[se_col]] <= true_values[p] & 
      true_values[p] <= res_df[[est_col]] + 1.645 * res_df[[se_col]], na.rm = TRUE)
    metrics$CP95[metrics$Parameter == p] <- mean(
      res_df[[est_col]] - 1.96 * res_df[[se_col]] <= true_values[p] & 
      true_values[p] <= res_df[[est_col]] + 1.96 * res_df[[se_col]], na.rm = TRUE)
  }
  
  cat("\n=== Monte Carlo Results (", n_sim, " replications) ===\n")
  print(metrics, digits = 5)
  
  invisible(list(metrics = metrics, raw_results = res_df))
}


# Run the simulation
mc <- monte_carlo_boxcox(n_sim = 1000, n = 2000,alpha_true=1.0, rho = 0.12)   

