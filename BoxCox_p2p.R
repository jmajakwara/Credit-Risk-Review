 
# ============================================================
# BOX-COX TRANSFORMATION CURE MODEL
# Lognormal
# ============================================================
rm(list=ls())
options(encoding = "UTF-8")
library(survival)
library(ggplot2)
library(MASS)
library(caret)
library(rsample)
library(dplyr)
library(data.table)
library(flexsurv)
library(flexsurvcure)
library(broom)
library(tidyverse)
library(censored)
library(parsnip)
library(tidymodels)

library(numDeriv)

# ===============================
# DATA PREPARATION 
# ===============================
loans  <- read.csv("loans_final_jj.csv",header=TRUE,sep=",",na.strings = "NA")

loans <- loans %>%
  dplyr::mutate(
    dplyr::across(c(time, OpenRevolvingMonthlyPayment, EmploymentStatusDuration,
             AmountDelinquent, RevolvingCreditBalance, AvailableBankcardCredit,
             LoanOriginalAmount), as.numeric),
    dplyr::across(c(ListingCategory, EmploymentStatus, CreditScoreRange,Term, Occupation, IsBorrowerHomeowner,CurrentlyInGroup), 
           as.factor)
  )

loans <- loans %>% dplyr::filter(DebtToIncomeRatio < 1 | is.na(DebtToIncomeRatio))

loans <- loans %>% dplyr::select(-2,-Occupation,-ProsperScore)

loans_cc <- loans %>% na.omit()

dim(loans)
dim(loans_cc)
#str(loans_cc)

loans_data <- loans_cc |> dplyr::select(time,status,DebtToIncomeRatio,Term,BorrowerRate,CreditScoreRange,
                                      BankcardUtilization,InquiriesLast6Months,OpenCreditLines,IsBorrowerHomeowner)

cat("Fitting standard survival models \n")

models <- list()

# Exponential
models$Exponential <- survival_reg(dist = "exp") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ BorrowerRate + OpenCreditLines +  InquiriesLast6Months + DebtToIncomeRatio +  BankcardUtilization + CreditScoreRange,
				 data = loans_data)

# Weibull
models$Weibull <- survival_reg(dist = "weibull") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines, data = loans_data)

# Gamma
models$Gamma <- survival_reg(dist = "gamma") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines+IsBorrowerHomeowner, data = loans_data)

# Lognormal
models$Lognormal <- survival_reg(dist = "lnorm") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines+IsBorrowerHomeowner, data = loans_data)

# Log-logistic
models$Loglogistic <- survival_reg(dist = "llogis") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines+IsBorrowerHomeowner, data = loans_data)

# Gompertz
models$Gompertz <- survival_reg(dist = "gompertz") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") |>
  fit(Surv(time, status) ~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines+IsBorrowerHomeowner, data = loans_data)

# Extract metrics
results <- lapply(names(models), function(name) {
  fit <- models[[name]]
  glance_fit <- glance(fit)
  data.frame(
    Model = name,
    LogLik = round(glance_fit$logLik, 2),
    AIC = round(glance_fit$AIC, 2),
    BIC = round(glance_fit$BIC, 2),
    n = glance_fit$N,
    df = glance_fit$df
  )
}) |> bind_rows()

print(results)




# ============================================================
# Negative Log-Likelihood
# ============================================================

boxcox_nll <- function(par, time, status, X) {
  
  alpha_raw <- par[1]
  alpha <- plogis(alpha_raw)
  
  beta <- par[2:(ncol(X) + 1)]
  eta  <- as.vector(X %*% beta)
  
  meanlog <- par[ncol(X) + 2]
  log_sdlog <- par[length(par)]
  sdlog   <- exp(log_sdlog)
  
  F_t <- plnorm(time, meanlog = meanlog, sdlog = sdlog)
  f_t <- dlnorm(time, meanlog = meanlog, sdlog = sdlog)
  
  if (alpha < 1e-5) {
    phi <- exp(eta)
  } else {
    phi <- exp(eta) / (1 + alpha * exp(eta))
  }
  
  if (alpha < 1e-5) {
    S_p <- exp(-phi * F_t)
    f_p <- phi * f_t * exp(-phi * F_t)
  } else {
    term <- 1 - alpha * phi * F_t
    term <- pmax(term, 1e-12)
    S_p  <- term^(1 / alpha)
    f_p  <- phi * f_t * term^(1 / alpha - 1)
  }
  
  S_p <- pmax(S_p, 1e-12)
  f_p <- pmax(f_p, 1e-12)
  
  ll <- sum(status * log(f_p) + (1 - status) * log(S_p))
  return(-ll)
}


# ============================================================
# Fit Function
# ============================================================

fit_boxcox <- function(data, maxit = 10000) {
  
  X <- model.matrix(~ DebtToIncomeRatio+Term+BorrowerRate+CreditScoreRange+BankcardUtilization+
                      InquiriesLast6Months+OpenCreditLines+IsBorrowerHomeowner, data = data)
  
  time   <- data$time
  status <- data$status
  
  start <- c(qlogis(0.4), rep(0, ncol(X)), 5.3, log(1.2))
  
  fit_optim <- optim(par = start, 
                     fn = boxcox_nll,
                     time = time, status = status, X = X,
                     method = "BFGS",
                     control = list(maxit = maxit, reltol = 1e-10),
                     hessian = TRUE)
  
  est <- fit_optim$par
  alpha_raw <- est[1]
  alpha_hat <- plogis(alpha_raw)
  
  names(est) <- c("alpha_raw", colnames(X), "meanlog", "log_sdlog")
  
 

H <- hessian(boxcox_nll,
             fit_optim$par,
             time=time,
             status=status,
             X=X)
  #V <- tryCatch(solve(H), error = function(e) NULL) 
  V <- tryCatch(solve(fit_optim$hessian), error = function(e) NULL) 

  se <- if(!is.null(V)) sqrt(diag(V)) else rep(NA, length(est))
  
  # Delta method SE for alpha = plogis(alpha_raw)
  se_alpha <- alpha_hat * (1 - alpha_hat) * se[1]
  
  # t-stats and p-values
  t_stat <- est / se
  p_value <- 2 * (1 - pt(abs(t_stat), df = nrow(data) - length(est)))
  p_value[is.na(se) | se < 1e-12] <- NA
  
  results <- data.frame(
    Parameter = names(est),
    Estimate  = est,
    Std.Error = se,
    Statistic = t_stat,
    Pvalue    = p_value
  )
  
  # Add transformed alpha row
  alpha_row <- data.frame(
    Parameter = "alpha",
    Estimate  = alpha_hat,
    Std.Error = se_alpha,
    Statistic = alpha_hat / se_alpha,
    Pvalue    = 2 * (1 - pt(abs(alpha_hat / se_alpha), df = nrow(data) - length(est)))
  )
  
  results <- rbind(alpha_row, results)
  rownames(results) <- NULL
  
  # Cure fraction
  eta_x <- X %*% est[2:(ncol(X)+1)]
  if (alpha_hat < 1e-5) {
    phi_x <- exp(eta_x)
    cure_prob <- exp(-phi_x)
  } else {
    phi_x <- exp(eta_x) / (1 + alpha_hat * exp(eta_x))
    cure_prob <- (1 - alpha_hat * phi_x)^(1 / alpha_hat)
  }
  
  list(
    results       = results,
    alpha         = alpha_hat,
    se_alpha      = se_alpha,
    cure_fraction = as.vector(cure_prob),
    logLik        = -fit_optim$value,
    AIC           = 2*length(est) - 2*(-fit_optim$value),
    BIC           = -2*(-fit_optim$value)+log(nrow(data))*length(est),
    convergence   = fit_optim$convergence
  )
}

# ============================================================
# Example
# ============================================================

fit <- fit_boxcox(loans_data)
fit$results
fit$AIC
fit$BIC
summary(fit$cure_fraction)
quantile(fit$cure_fraction, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
boxplot(fit$cure_fraction,horizontal = TRUE)


# Fitting mixture cure rate model
log_fit <- flexsurvcure(Surv(time_to_event,status) ~ CreditScoreRange,
                           anc=list(meanlog = ~ DebtToIncomeRatio + Term + BorrowerRate + BankcardUtilization +
                      InquiriesLast6Months + OpenCreditLines + IsBorrowerHomeowner),  
			 data = loans_data,
			 link = "logistic",   # incidence (cure probability) on logit scale
  			 dist = "lnorm",     # latency = log-logistic AFT
  			 mixture = TRUE       
)

print(tidy(log_fit))
print(glance(log_fit))

