###############################################################################
# Monte Carlo simulation for the lognormal Box-Cox transformation cure model
###############################################################################
rm(list = ls())
options(encoding = "UTF-8")
library(dplyr)
library(numDeriv)

# ---------------- OPTIONS ----------------
n_sim     <- 1000
n_obs     <- 1000          # 1000 or 2000
alpha_set <- 0.5           # true alpha for the single run below
rho       <- 0.12
run_all   <- TRUE         # TRUE: loop over alpha in {0, .25, .5, .75, 1} and n in {1000, 2000} and write CSV files

# True parameters. Order: intercept, X1B, X1C, X1D, X2, X3, X4, X5.

BETA_TRUE   <- c(`(Intercept)` = 1.6, X1B = -0.7, X1C = -0.49, X1D = -0.41, X2 = -5.65, X3 = -0.55, X4 = -0.06, X5 = 0.4)
MEANLOG_TRUE <- 3.2
SDLOG_TRUE   <- 1.3

# ---------------- DATA GENERATION ----------------

generate_logn_boxcox_data <- function(n = 1000, seed = 2000, rho = 0.12, alpha_true = 0.5, cens_rate = 0.001) {
  set.seed(seed)
  cor_matrix <- matrix(rho, 5, 5); diag(cor_matrix) <- 1
  Y <- MASS::mvrnorm(n = n, mu = rep(0, 5), Sigma = cor_matrix)        
  U <- as.data.frame(pnorm(Y)); colnames(U) <- paste0("U", 1:5)

  X <- U |>
    dplyr::mutate(
      X1 = case_when(U1 <= 0.04 ~ "A", U1 <= 0.24 ~ "B", U1 <= 0.86 ~ "C", TRUE ~ "D"),   # 4%, 20%, 62%, 14%
      X2 = pmin(pmax(U2 * 0.36, 0.04), 0.36),                                 # uniform(0, .36) clamped to [.04, .36]
      X3 = pmax(pmin(qlnorm(U3, meanlog = -1.555, sdlog = 0.58), 0.99), 0.01),
      X4 = pmin(rpois(n, lambda = 15 * U4), 22),                              # Poisson with random rate 15*U4, capped at 22
      X5 = pmin(pmax(U5 * 2.5, 0), 2.5)
    ) |>
    dplyr::mutate(X1 = factor(X1, levels = c("A", "B", "C", "D")))
  X_mat <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = X)
  stopifnot(all(colnames(X_mat) == names(BETA_TRUE)))

  eta <- as.vector(X_mat %*% BETA_TRUE)
  if (alpha_true < 1e-5) { phi <- exp(eta); cure_prob <- exp(-phi)
  } else { phi <- exp(eta) / (1 + alpha_true * exp(eta)); cure_prob <- (1 - alpha_true * phi)^(1 / alpha_true) }


  C <- rexp(n, rate = cens_rate)
  time_obs <- numeric(n); status <- numeric(n); cured_true <- numeric(n)
  for (i in 1:n) {
    p0i <- cure_prob[i]
    if (runif(1) <= p0i) { time_obs[i] <- C[i]; status[i] <- 0; cured_true[i] <- 1
    } else {
      u2 <- runif(1)
      arg <- if (alpha_true < 1e-5) -log(p0i + (1 - p0i) * u2) / phi[i]
             else (1 - (p0i + (1 - p0i) * u2)^alpha_true) / (alpha_true * phi[i])   # F(t_e) from S_p = p0 + (1-p0) u2
      t_event <- exp(MEANLOG_TRUE + SDLOG_TRUE * qnorm(arg))
      time_obs[i] <- min(t_event, C[i]); status[i] <- as.numeric(t_event <= C[i])
    }
  }
  data.frame(time = time_obs, status = status, X1 = X$X1, X2 = X$X2, X3 = X$X3, X4 = X$X4, X5 = X$X5, cured_true = cured_true)
}

# ---------------- LIKELIHOOD (unchanged) ----------------

boxcox_nll <- function(par, time, status, X) {
  alpha <- plogis(par[1])
  beta <- par[2:(ncol(X) + 1)]
  eta <- as.vector(X %*% beta)
  meanlog <- par[ncol(X) + 2]
  sdlog <- exp(par[length(par)])
  F_t <- plnorm(time, meanlog = meanlog, sdlog = sdlog)
  f_t <- dlnorm(time, meanlog = meanlog, sdlog = sdlog)
  if (alpha < 1e-5) {
    phi <- exp(eta); S_p <- exp(-phi * F_t); f_p <- phi * f_t * exp(-phi * F_t)
  } else {
    phi <- exp(eta) / (1 + alpha * exp(eta))
    term <- pmax(1 - alpha * phi * F_t, 1e-12)
    S_p <- term^(1 / alpha); f_p <- phi * f_t * term^(1 / alpha - 1)
  }
  S_p <- pmax(S_p, 1e-12); f_p <- pmax(f_p, 1e-12)
  -sum(status * log(f_p) + (1 - status) * log(S_p))
}

# ---------------- MONTE CARLO ----------------

monte_carlo_boxcox <- function(n_sim = 1000, n = 1000, rho = 0.12, alpha_true = 0.5, verbose = TRUE) {
  results_list <- vector("list", n_sim); ok <- logical(n_sim); se_ok <- logical(n_sim)
  for (i in 1:n_sim) {
    if (verbose) cat("Simulation", i, "of", n_sim, "\r")
    dat <- generate_logn_boxcox_data(n = n, seed = 1000 + i, alpha_true = alpha_true, rho = rho)
    X <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = dat)
    start <- c(qlogis(0.5), rep(0, ncol(X)), 3.5, log(1.25))                 # FIX: log scale for sdlog
    fit <- tryCatch(optim(start, boxcox_nll, time = dat$time, status = dat$status, X = X,
                          method = "BFGS", control = list(maxit = 10000, reltol = 1e-10)), error = function(e) NULL)
    k <- ncol(X)
    if (is.null(fit) || fit$convergence != 0 || !is.finite(fit$value)) { results_list[[i]] <- rep(NA_real_, 2 * (k + 3)); next }
    ok[i] <- TRUE
    est <- fit$par; alpha_hat <- plogis(est[1])
    H <- tryCatch(numDeriv::hessian(boxcox_nll, est, time = dat$time, status = dat$status, X = X), error = function(e) NULL)
    V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
    if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) > 0)) {
      se <- sqrt(diag(V)); se_ok[i] <- TRUE
    } else se <- rep(NA_real_, length(est))
    results_list[[i]] <- c(alpha_hat, est[-1], alpha_hat * (1 - alpha_hat) * se[1], se[-1])
  }
  param_names <- c("alpha", colnames(X), "meanlog", "log_sdlog")
  res_mat <- do.call(rbind, results_list)
  colnames(res_mat) <- c(paste0(param_names, "_est"), paste0(param_names, "_se"))
  res_df <- as.data.frame(res_mat)

  true_values <- c(alpha = alpha_true, BETA_TRUE, meanlog = MEANLOG_TRUE, log_sdlog = log(SDLOG_TRUE))   
  tv <- true_values[param_names]
  est_mat <- as.matrix(res_df[, paste0(param_names, "_est")])
  err <- sweep(est_mat, 2, tv)                                                
  mean_est <- colMeans(est_mat, na.rm = TRUE)
  metrics <- data.frame(Parameter = param_names, True = unname(tv), MeanEst = unname(mean_est),
                        Bias = unname(mean_est - tv),
                        PctBias = ifelse(tv == 0, NA_real_, 100 * unname(mean_est - tv) / abs(unname(tv))),
                        MSE = unname(colMeans(err^2, na.rm = TRUE)), RMSE = sqrt(unname(colMeans(err^2, na.rm = TRUE))))
  cover <- function(z) sapply(param_names, function(p) {
    e <- res_df[[paste0(p, "_est")]]; s <- res_df[[paste0(p, "_se")]]
    mean(e - z * s <= tv[p] & tv[p] <= e + z * s, na.rm = TRUE) })
  metrics$CP90 <- unname(cover(1.645)); metrics$CP95 <- unname(cover(1.96))

  cat(sprintf("\n=== Monte Carlo results: alpha = %.2f, n = %d, %d replications (%d converged, %d with valid SE) ===\n",
              alpha_true, n, n_sim, sum(ok), sum(se_ok)))
  print(metrics, digits = 5, row.names = FALSE)
  invisible(list(metrics = metrics, raw_results = res_df, n_converged = sum(ok), n_valid_se = sum(se_ok)))
}

# ---------------- RUN ----------------

if (!run_all) {
  mc <- monte_carlo_boxcox(n_sim = n_sim, n = n_obs, rho = rho, alpha_true = alpha_set)
} else {
  for (nn in c(1000, 2000)) for (aa in c(0, 0.25, 0.5, 0.75, 1)) {
    mc <- monte_carlo_boxcox(n_sim = n_sim, n = nn, rho = rho, alpha_true = aa, verbose = FALSE)
    write.csv(mc$metrics, sprintf("mc_boxcox_alpha%.2f_n%d.csv", aa, nn), row.names = FALSE)
    saveRDS(mc$raw_results, sprintf("mc_boxcox_alpha%.2f_n%d_raw.rds", aa, nn))
  }
}