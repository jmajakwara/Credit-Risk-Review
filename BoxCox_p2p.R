# ============================================================
# MIXTURE CURE MODEL (flexsurvcure) AND BOX-COX TRANSFORMATION (BCT) CURE MODEL, LOGNORMAL LATENCY
# cleaned Prosper data; SurvSHAP(t) for the mixture cure model only
#
#   Mixture cure (logistic link):  S_p(t|x,z) = pi(x) + {1 - pi(x)} S_u(t|z),   logit pi(x) = x'b      (pi = cure probability)
#   BCT cure model:                S_p(t|x,z) = {1 - alpha*phi*F(t|z)}^(1/alpha),  phi = exp(eta)/(1 + alpha*exp(eta)),
#                                  alpha = plogis(alpha_raw) in (0,1);  alpha = 1 is the logistic mixture cure model
#                                  (then eta = -logit pi, i.e. b_BCT = -b_mixture); alpha -> 0 is the promotion-time model.
#   Latency (both): lognormal AFT, F(t|z) = F0(t exp(-z'gamma)),  gamma > 0 = later default (time ratio exp(gamma)).
#
# SIGNS: mixture cure coefficients are on the log-odds of CURE (positive = MORE cure);
#        BCT cure coefficients are on eta (positive = LESS cure).
#
# OPTIONS (edit the block below)
#   term_keep : 60   -> 60-month loans only (the only term with enough follow-up for a Maller-Zhou plateau)
#               NULL -> all loan terms (Term dummies are added as covariates)
#   design    : "split"   -> CURE: CreditScoreRange only;  LATENCY (AFT): all other covariates
#               "allcure" -> every covariate in the cure part, latency has no covariates
#   run_standard_models : fit no-cure parametric survival models (flexsurv) for comparison
#   run_shap  : SurvSHAP(t) for the MIXTURE cure model (slow)
#
# Estimates, SE, z and p-values are reported on the ORIGINAL covariate scale for both models
# (fitting is done on standardised covariates for numerical stability, then back-transformed with the delta method).
# The cure intercept refers to all cure covariates = 0 (original scale); the baseline meanlog refers to the
# latency covariates at their MEAN.
#
# Data: prosper_analysis.csv (cleaned file).  default = 1 only for event_type 1; matured/active (2) and early repaid (3) are censored.
# ============================================================

rm(list = ls())
options(encoding = "UTF-8")
library(dplyr)
library(survival)
library(numDeriv)
library(ggplot2)
library(flexsurvcure)

# ---------------- OPTIONS ----------------
data_file           <- "prosper_analysis.csv"          # <- path to the cleaned file (all loan terms; term_keep selects the term)
term_keep           <- 60                              # 60 or NULL (all terms)
design              <- "split"                         # "split" or "allcure"
dti_max             <- 1                               # keep DebtToIncomeRatio < dti_max (Inf keeps all)
run_standard_models <- TRUE
run_shap            <- TRUE                            # SurvSHAP(t) for the mixture cure model
shap_n              <- 50                              # test-set loans explained by SurvSHAP(t); Inf = all test loans (very slow)
shap_times          <- "grid"                          # "grid": 100 equally spaced times (fast);  "unique": every distinct test time (slow)
out_prefix          <- paste0("mcm_bct_lognormal_", if (is.null(term_keep)) "allterms" else paste0("term", term_keep), "_", design)
stopifnot(design %in% c("split", "allcure"))

# ===============================
# DATA PREPARATION
# ===============================
loans <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
need <- c("time", "BorrowerRate", "CreditScoreRange", "DebtToIncomeRatio", "Term",
          "BankcardUtilization", "InquiriesLast6Months", "OpenCreditLines", "IsBorrowerHomeowner")
stopifnot(all(need %in% names(loans)))
if ("event_type" %in% names(loans)) {
  cat("Event types in the file (1 default, 2 matured/active, 3 early repayment):\n"); print(table(loans$event_type))
  loans$status <- as.integer(loans$event_type == 1L)
  if ("default" %in% names(loans)) stopifnot(all(loans$status == loans$default))
} else if ("default" %in% names(loans)) loans$status <- as.integer(loans$default) else loans$status <- as.integer(as.character(loans$status))

loans_data <- loans %>%
  dplyr::filter(is.na(DebtToIncomeRatio) | DebtToIncomeRatio < dti_max) %>%
  dplyr::select(dplyr::all_of(c("time", "status", setdiff(need, "time")))) %>%
  dplyr::mutate(
    time = as.numeric(time),
    Term = factor(Term),
    CreditScoreRange = factor(as.character(CreditScoreRange), levels = c("Excellent", "Very good", "Good", "Fair")),  # Excellent = reference
    IsBorrowerHomeowner = factor(as.character(IsBorrowerHomeowner), levels = c("FALSE", "TRUE"))
  ) %>%
  stats::na.omit() %>%
  dplyr::filter(time > 0) %>%
  { if (is.null(term_keep)) . else dplyr::filter(., Term == as.character(term_keep)) } %>%
  droplevels()
stopifnot(all(loans_data$status %in% 0:1), nlevels(loans_data$CreditScoreRange) >= 2)
cat(sprintf("\nn = %d | defaults = %d (%.1f%%) | censored = %d | max follow-up = %.0f days | terms: %s\n",
            nrow(loans_data), sum(loans_data$status), 100 * mean(loans_data$status),
            sum(loans_data$status == 0), max(loans_data$time), paste(levels(loans_data$Term), collapse = ", ")))

# covariate formulas for the two parts of the model
term_part  <- if (nlevels(loans_data$Term) > 1) "Term + " else ""
all_cov    <- paste0("DebtToIncomeRatio + ", term_part, "BorrowerRate + CreditScoreRange + BankcardUtilization + ",
                     "InquiriesLast6Months + OpenCreditLines + IsBorrowerHomeowner")
cure_formula <- if (design == "split") ~ CreditScoreRange else as.formula(paste("~", all_cov))
lat_formula  <- if (design == "split")
  as.formula(paste0("~ DebtToIncomeRatio + ", term_part, "BorrowerRate + BankcardUtilization + ",
                    "InquiriesLast6Months + OpenCreditLines + IsBorrowerHomeowner")) else ~ 1

# ============================================================
# NO-CURE BENCHMARKS (flexsurv). Each distribution keeps the covariate set of your original script.
#   AIC / BIC compare models fitted to the same data by maximum likelihood, so the covariate sets need not be equal;
#   only a likelihood-ratio test would require nested models.
# ============================================================
# Weibull: flexsurv takes its starting values from survreg(), which does not converge on raw covariates
# -> non-finite initial objective.  Own ML on standardised covariates supplies the starting values.
weibull_start <- function(ff, data) {
  mf <- model.frame(ff, data); y <- model.response(mf); X <- model.matrix(ff, mf)
  tt <- y[, 1]; ev <- y[, 2] == 1
  mu <- colMeans(X[, -1, drop = FALSE]); sdv <- apply(X[, -1, drop = FALSE], 2, sd)
  Xs <- X; Xs[, -1] <- sweep(sweep(X[, -1, drop = FALSE], 2, mu), 2, sdv, "/")
  lt <- log(tt)
  nll <- function(p) { k <- exp(p[1]); ls <- as.vector(Xs %*% p[-1]); z <- k * (lt - ls)
    v <- -(sum((log(k) - ls + (k - 1) * (lt - ls) - exp(z))[ev]) - sum(exp(z)[!ev])); if (is.finite(v)) v else 1e10 }
  sh0 <- 1.64 / var(lt)
  o <- optim(c(log(sh0), mean(lt) + 0.572 / sh0, rep(0, ncol(X) - 1)), nll, method = "BFGS",
             control = list(maxit = 5000, reltol = 1e-14))
  b <- o$par[-1]; b_orig <- c(b[1] - sum(b[-1] * mu / sdv), b[-1] / sdv)
  list(inits = c(exp(o$par[1]), exp(b_orig[1]), b_orig[-1]), loglik = -o$value, npar = length(o$par), code = o$convergence)
}

results_std <- NULL
if (run_standard_models) {
  library(flexsurv)
  tp <- if (nlevels(loans_data$Term) > 1) "Term + " else ""
  rhs <- list(
    exp      = "BorrowerRate + OpenCreditLines + InquiriesLast6Months + DebtToIncomeRatio + BankcardUtilization + CreditScoreRange",
    weibull  = paste0("DebtToIncomeRatio + ", tp, "BorrowerRate + CreditScoreRange + BankcardUtilization + InquiriesLast6Months + OpenCreditLines"),
    gamma    = all_cov, lnorm = all_cov, llogis = all_cov, gompertz = all_cov, gengamma = all_cov)
  std <- lapply(names(rhs), function(d) {
    ff <- as.formula(paste("Surv(time, status) ~", rhs[[d]]))
    w <- if (d == "weibull") weibull_start(ff, loans_data) else NULL
    tryCatch({
      m <- if (is.null(w)) flexsurv::flexsurvreg(ff, data = loans_data, dist = d)       # inits = NULL is not accepted by flexsurvreg
           else flexsurv::flexsurvreg(ff, data = loans_data, dist = d, inits = w$inits)
      data.frame(Model = d, LogLik = m$loglik, npar = m$npars, AIC = -2 * m$loglik + 2 * m$npars,
                 BIC = -2 * m$loglik + log(m$N) * m$npars)
    }, error = function(e) {
      message("flexsurv ", d, " failed: ", conditionMessage(e))
      if (!is.null(w))                                      # direct ML result for the Weibull
        data.frame(Model = paste0(d, " (direct ML)"), LogLik = w$loglik, npar = w$npar,
                   AIC = -2 * w$loglik + 2 * w$npar, BIC = -2 * w$loglik + log(nrow(loans_data)) * w$npar) else NULL
    })
  })
  results_std <- do.call(rbind, std)
  cat("\n========== NO-CURE PARAMETRIC MODELS (flexsurv; original covariate sets, same data) ==========\n")
  print(results_std, digits = 6, row.names = FALSE)
}

# ============================================================
# DESIGN MATRICES (stored so that prediction on new data uses the same coding and scaling)
# ============================================================
build_design <- function(data, cure_f, lat_f, xlev, sc = NULL) {
  mf_x <- model.frame(cure_f, data = data, xlev = xlev[intersect(names(xlev), all.vars(cure_f))], na.action = na.pass)
  X <- model.matrix(cure_f, mf_x)
  Z <- if (length(all.vars(lat_f)) > 0) {
    mf_z <- model.frame(lat_f, data = data, xlev = xlev[intersect(names(xlev), all.vars(lat_f))], na.action = na.pass)
    model.matrix(lat_f, mf_z)[, -1, drop = FALSE]
  } else matrix(0, nrow(data), 0)
  if (is.null(sc)) {
    sc <- list(mx = colMeans(X[, -1, drop = FALSE]), sx = apply(X[, -1, drop = FALSE], 2, sd),
               mz = colMeans(Z), sz = if (ncol(Z)) apply(Z, 2, sd) else numeric(0))
    if (any(sc$sx == 0) || any(sc$sz == 0)) stop("constant covariate column")
  }
  Xs <- X; Xs[, -1] <- sweep(sweep(X[, -1, drop = FALSE], 2, sc$mx), 2, sc$sx, "/")
  Zs <- if (ncol(Z)) sweep(sweep(Z, 2, sc$mz), 2, sc$sz, "/") else Z
  list(X = X, Z = Z, Xs = Xs, Zs = Zs, sc = sc)
}

# back-transform  e = (beta_1..beta_k, meanlog, log sdlog, gamma_1..gamma_q)  from the standardised to the original scale
#   beta_orig = A beta_std,  gamma_orig = gamma_std / sd,  sdlog = exp(log sdlog);  SE by the delta method
back_transform <- function(e, V, D) {
  k <- ncol(D$Xs); q <- ncol(D$Zs); sc <- D$sc
  A <- diag(1 / c(1, sc$sx), k); A[1, -1] <- -sc$mx / sc$sx
  sdl <- exp(e[k + 2])
  Jm <- diag(length(e)); Jm[seq_len(k), seq_len(k)] <- A; Jm[k + 2, k + 2] <- sdl
  if (q) Jm[k + 2 + seq_len(q), k + 2 + seq_len(q)] <- diag(1 / sc$sz, q)
  est <- c(as.vector(A %*% e[seq_len(k)]), e[k + 1], sdl, if (q) e[k + 2 + seq_len(q)] / sc$sz)
  se  <- if (!is.null(V)) sqrt(pmax(diag(Jm %*% V %*% t(Jm)), 0)) else rep(NA_real_, length(est))
  list(est = est, se = se)
}
param_names <- function(D) c(colnames(D$X), "meanlog", "sdlog", if (ncol(D$Z)) paste0("time-ratio: ", colnames(D$Z)))

safe_vcov <- function(H) {
  min_eig <- min(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
  V <- if (is.finite(min_eig) && min_eig > 0) tryCatch(solve(H), error = function(e) NULL) else NULL
  list(V = V, min_eig = min_eig)
}

# ============================================================
# (A) MIXTURE CURE MODEL
#   direct ML (own likelihood) -> starting values for flexsurvcure::flexsurvcure(); flexsurvcure gives the reported fit.
#   If flexsurvcure fails, the direct ML fit is reported instead (flagged in the output).
# ============================================================
mixture_nll <- function(par, time, status, Xs, Zs) {
  k <- ncol(Xs); q <- ncol(Zs)
  lin <- pmin(pmax(as.vector(Xs %*% par[seq_len(k)]), -40), 40)       # logit of the cure probability
  ml <- par[k + 1]; sdl <- exp(par[k + 2])
  s  <- if (q) as.vector(Zs %*% par[k + 2 + seq_len(q)]) else 0
  ye <- time * exp(-s)
  lf <- dlnorm(ye, ml, sdl, log = TRUE) - s
  la <- plogis(lin, log.p = TRUE)                                      # log pi
  lb <- plogis(lin, lower.tail = FALSE, log.p = TRUE)                  # log (1 - pi)
  lS <- plnorm(ye, ml, sdl, lower.tail = FALSE, log.p = TRUE)
  b  <- lb + lS; m <- pmax(la, b); lcens <- m + log(exp(la - m) + exp(b - m))
  ev <- status == 1
  ll <- sum((lb + lf)[ev]) + sum(lcens[!ev])
  if (is.finite(ll)) -ll else 1e10
}

fit_mixture <- function(data, cure_f = cure_formula, lat_f = lat_formula, maxit = 3000, verbose = TRUE) {
  xlev <- lapply(data[sapply(data, is.factor)], levels)
  D <- build_design(data, cure_f, lat_f, xlev)
  Xs <- D$Xs; Zs <- D$Zs; k <- ncol(Xs); q <- ncol(Zs)
  time <- data$time; status <- data$status; n <- length(time)

  # step 1: direct ML
  p0km <- min(max(min(survfit(Surv(time, status) ~ 1)$surv), 0.05), 0.95)
  lt <- log(time[status == 1])
  st <- c(qlogis(p0km), rep(0, k - 1), mean(lt), log(sd(lt)), rep(0, q))
  o <- optim(st, mixture_nll, time = time, status = status, Xs = Xs, Zs = Zs, method = "BFGS",
             control = list(maxit = maxit, reltol = 1e-12))
  if (verbose) cat(sprintf("  direct ML (mixture): logLik %.3f (code %d)\n", -o$value, o$convergence))

  # step 2: flexsurvcure on the standardised covariates (cure covariates -> theta; latency covariates -> meanlog via anc)
  cn <- paste0("c_", make.names(colnames(Xs)[-1])); zn <- if (q) paste0("z_", make.names(colnames(Zs))) else character(0)
  dat <- data.frame(time = time, status = status)
  dat <- cbind(dat, setNames(as.data.frame(Xs[, -1, drop = FALSE]), cn))
  if (q) dat <- cbind(dat, setNames(as.data.frame(Zs), zn))
  f   <- as.formula(paste("Surv(time, status) ~", paste(cn, collapse = " + ")))
  anc <- if (q) list(meanlog = as.formula(paste("~", paste(zn, collapse = " + ")))) else NULL
  inits <- c(plogis(o$par[1]), o$par[k + 1], exp(o$par[k + 2]), o$par[1 + seq_len(k - 1)], o$par[k + 2 + seq_len(q)])
  nm <- c("theta", cn, "meanlog", "sdlog", if (q) paste0("meanlog(", zn, ")"))
  m <- tryCatch(flexsurvcure::flexsurvcure(f, data = dat, dist = "lnorm", link = "logistic", mixture = TRUE,
                                           anc = anc, inits = inits, control = list(maxit = 5000, reltol = 1e-12)),
                error = function(e) { message("flexsurvcure failed: ", conditionMessage(e)); NULL })
  source <- "flexsurvcure"
  if (!is.null(m) && all(nm %in% rownames(m$res.t)) && is.finite(m$loglik)) {
    # flexsurvcure parameter order: theta (intercept, logit scale), cure slopes, meanlog, log sdlog, latency effects
    e  <- m$res.t[nm, "est"]; V <- vcov(m)[nm, nm, drop = FALSE]
    # reorder to (beta_1..beta_k, meanlog, log sdlog, gamma)
    e  <- c(e["theta"], e[cn], e["meanlog"], e["sdlog"], if (q) e[paste0("meanlog(", zn, ")")])
    ord <- nm[c(1, seq_along(cn) + 1, length(cn) + 2, length(cn) + 3, if (q) length(cn) + 3 + seq_len(q))]
    V  <- V[ord, ord, drop = FALSE]
    ll <- m$loglik; npar <- m$npars; min_eig <- min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
    e <- unname(e); V <- unname(V)
    if (verbose) cat(sprintf("  flexsurvcure: logLik %.3f\n", ll))
  } else {
    if (!is.null(m)) message("flexsurvcure output names not as expected; reporting the direct ML fit")
    source <- "direct ML (flexsurvcure unavailable)"
    e <- o$par; ll <- -o$value; npar <- length(e)
    H <- numDeriv::hessian(mixture_nll, e, time = time, status = status, Xs = Xs, Zs = Zs)
    sv <- safe_vcov(H); V <- sv$V; min_eig <- sv$min_eig
    if (is.null(V)) warning("Hessian not positive definite: standard errors unavailable")
  }
  bt <- back_transform(e, V, D)
  tab <- data.frame(Parameter = param_names(D), Estimate = bt$est, Std.Error = bt$se)
  tab$z <- tab$Estimate / tab$Std.Error; tab$Pvalue <- 2 * pnorm(-abs(tab$z)); rownames(tab) <- NULL
  cure_i <- plogis(as.vector(Xs %*% e[seq_len(k)]))
  model <- list(type = "mixture", source = source, par = e, k = k, q = q,
                meanlog = e[k + 1], sdlog = exp(e[k + 2]),
                beta = e[seq_len(k)], gamma = if (q) e[k + 2 + seq_len(q)] else numeric(0),
                cure_formula = cure_f, lat_formula = lat_f, xlev = xlev, sc = D$sc,
                results = tab, cure_fraction = cure_i, logLik = ll, AIC = 2 * npar - 2 * ll,
                BIC = log(n) * npar - 2 * ll, n_par = npar, min_hessian_eigen = min_eig,
                flexsurvcure_object = m)
  class(model) <- c("mixture_cure_model", "cure_fit")
  model
}

# ============================================================
# (B) BOX-COX TRANSFORMATION CURE MODEL
#   par = (alpha_raw, beta_1..beta_k, meanlog, log sdlog, gamma_1..gamma_q)  on standardised covariates
#   events:   log f_p = log(phi) + log f(t_e) + (1/alpha - 1) log(1 - alpha*phi*F(t_e)),  t_e = t exp(-z'gamma)
#   censored: log S_p = log(1 - alpha*phi*F(t_e)) / alpha
# ============================================================
boxcox_nll <- function(par, time, status, Xs, Zs) {
  k <- ncol(Xs); q <- ncol(Zs)
  a   <- plogis(par[1])
  eta <- pmin(pmax(as.vector(Xs %*% par[1 + seq_len(k)]), -30), 30)
  ml  <- par[k + 2]; sdl <- exp(par[k + 3])
  s   <- if (q) as.vector(Zs %*% par[k + 3 + seq_len(q)]) else 0
  ye  <- time * exp(-s)
  Ft  <- plnorm(ye, meanlog = ml, sdlog = sdl)
  lf  <- dlnorm(ye, meanlog = ml, sdlog = sdl, log = TRUE) - s
  if (a < 1e-5) {
    logfp <- eta + lf - exp(eta) * Ft
    logS  <- -exp(eta) * Ft
  } else {
    phi <- 1 / (a + exp(-eta))
    B   <- pmax(1 - a * phi * Ft, 1e-300)
    logfp <- log(phi) + lf + (1 / a - 1) * log(B)
    logS  <- log(B) / a
  }
  ev <- status == 1
  ll <- sum(logfp[ev]) + sum(logS[!ev])
  if (is.finite(ll)) -ll else 1e10
}

cure_prob_fun <- function(eta, a) if (a < 1e-5) exp(-exp(eta)) else exp(-log1p(a * exp(eta)) / a)

fit_boxcox <- function(data, cure_f = cure_formula, lat_f = lat_formula, maxit = 3000, verbose = TRUE) {
  xlev <- lapply(data[sapply(data, is.factor)], levels)
  D <- build_design(data, cure_f, lat_f, xlev)
  Xs <- D$Xs; Zs <- D$Zs; k <- ncol(Xs); q <- ncol(Zs)
  time <- data$time; status <- data$status; n <- length(time)

  p0km <- min(max(min(survfit(Surv(time, status) ~ 1)$surv), 0.05), 0.95)
  lt <- log(time[status == 1])
  best <- NULL
  for (a0 in c(0.1, 0.5, 0.9)) {
    b0 <- log((p0km^(-a0) - 1) / a0)
    st <- c(qlogis(a0), b0, rep(0, k - 1), mean(lt), log(sd(lt)), rep(0, q))
    o <- tryCatch(optim(st, boxcox_nll, time = time, status = status, Xs = Xs, Zs = Zs, method = "BFGS",
                        control = list(maxit = maxit, reltol = 1e-12)), error = function(e) NULL)
    if (verbose && !is.null(o)) cat(sprintf("  start alpha %.1f: logLik %.3f (alpha-hat %.3f, code %d)\n",
                                            a0, -o$value, plogis(o$par[1]), o$convergence))
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) stop("optimisation failed from all starting values")
  est <- best$par; a_hat <- plogis(est[1])

  H <- numDeriv::hessian(boxcox_nll, est, time = time, status = status, Xs = Xs, Zs = Zs)
  sv <- safe_vcov(H); V <- sv$V; min_eig <- sv$min_eig
  if (is.null(V)) warning("Hessian not positive definite: standard errors unavailable (alpha at a boundary?)")
  se_raw <- if (!is.null(V)) sqrt(diag(V)) else rep(NA_real_, length(est))

  bt <- back_transform(est[-1], if (!is.null(V)) V[-1, -1, drop = FALSE] else NULL, D)
  se_alpha <- if (!is.null(V)) a_hat * (1 - a_hat) * se_raw[1] else NA_real_          # delta method
  tab <- data.frame(Parameter = c("alpha", param_names(D)), Estimate = c(a_hat, bt$est), Std.Error = c(se_alpha, bt$se))
  tab$z <- tab$Estimate / tab$Std.Error; tab$Pvalue <- 2 * pnorm(-abs(tab$z)); rownames(tab) <- NULL

  eta_i <- as.vector(Xs %*% est[1 + seq_len(k)])
  cure_i <- cure_prob_fun(eta_i, a_hat)
  ll <- -best$value; npar <- length(est)
  model <- list(type = "boxcox", source = "direct ML", alpha = a_hat, par = est, k = k, q = q,
                meanlog = est[k + 2], sdlog = exp(est[k + 3]),
                beta = est[1 + seq_len(k)], gamma = if (q) est[k + 3 + seq_len(q)] else numeric(0),
                cure_formula = cure_f, lat_formula = lat_f, xlev = xlev, sc = D$sc,
                results = tab, cure_fraction = cure_i, logLik = ll, AIC = 2 * npar - 2 * ll,
                BIC = log(n) * npar - 2 * ll, n_par = npar, convergence = best$convergence,
                min_hessian_eigen = min_eig)
  class(model) <- c("boxcox_cure_model", "cure_fit")
  model
}

print.cure_fit <- function(x, ...) {
  if (x$type == "mixture")
    cat("\nMIXTURE CURE MODEL (logistic link; fitted by ", x$source, ")\n",
        "CURE part: coefficients on log-odds of CURE (positive = MORE cure).  ",
        if (x$q) "LATENCY: time ratios, positive = LATER default.\n" else "\n", sep = "")
  else
    cat("\nBOX-COX TRANSFORMATION CURE MODEL\nCURE part (eta): positive coefficient = LOWER cure probability.  ",
        if (x$q) "LATENCY: time ratios, positive = LATER default.\n" else "\n", sep = "")
  print(x$results, digits = 4, row.names = FALSE)
  cat(sprintf("\nlogLik %.2f | AIC %.2f | BIC %.2f | npar %d | min Hessian/vcov eigenvalue %.3g\n",
              x$logLik, x$AIC, x$BIC, x$n_par, x$min_hessian_eigen))
  if (x$type == "boxcox")
    cat(sprintf("alpha-hat = %.4f%s\n", x$alpha,
                if (x$alpha > 0.99 || x$alpha < 0.01) "  (at the boundary of (0,1): the Wald SE for alpha is not meaningful)" else ""))
  cat(sprintf("mean cure fraction = %.3f\n", mean(x$cure_fraction)))
  invisible(x)
}

# ============================================================
# PREDICTION (population survival S_p(t | x)) -- used for the hold-out check and by survex
# ============================================================
predict_boxcox_survival <- function(model, newdata, times) {
  D <- build_design(as.data.frame(newdata), model$cure_formula, model$lat_formula, model$xlev, model$sc)
  a <- model$alpha
  eta <- pmin(pmax(as.vector(D$Xs %*% model$beta), -30), 30)
  s <- if (model$q) as.vector(D$Zs %*% model$gamma) else rep(0, nrow(D$Xs))
  phi <- if (a < 1e-5) exp(eta) else 1 / (a + exp(-eta))
  S <- matrix(NA_real_, nrow(D$Xs), length(times))
  for (j in seq_along(times)) {
    Ft <- plnorm(times[j] * exp(-s), meanlog = model$meanlog, sdlog = model$sdlog)
    S[, j] <- if (a < 1e-5) exp(-phi * Ft) else pmax(1 - a * phi * Ft, 1e-300)^(1 / a)
  }
  S
}
predict_mixture_survival <- function(model, newdata, times) {
  D <- build_design(as.data.frame(newdata), model$cure_formula, model$lat_formula, model$xlev, model$sc)
  pc <- plogis(as.vector(D$Xs %*% model$beta))
  s  <- if (model$q) as.vector(D$Zs %*% model$gamma) else rep(0, nrow(D$Xs))
  S <- matrix(NA_real_, nrow(D$Xs), length(times))
  for (j in seq_along(times))
    S[, j] <- pc + (1 - pc) * plnorm(times[j] * exp(-s), meanlog = model$meanlog, sdlog = model$sdlog, lower.tail = FALSE)
  S
}

# ============================================================
# KAPLAN-MEIER CURVES STRATIFIED BY CREDIT SCORE (log-rank p-value and number-at-risk table; PDF)
# ============================================================
if (requireNamespace("survminer", quietly = TRUE)) {
  km_fit <- survfit(Surv(time, status) ~ CreditScoreRange, data = loans_data)
  km_plot <- survminer::ggsurvplot(km_fit, data = loans_data, pval = TRUE, risk.table = TRUE, censor = TRUE,
                                   xlab = "Time (Days)", ylab = "Survival Probability", legend.title = "Strata",
                                   title = "Kaplan-Meier Survival Curves stratified by Credit Score",
                                   ylim = c(0, 1), break.time.by = 400, ggtheme = theme_minimal())
  pdf(paste0(out_prefix, "_KM_by_creditscore.pdf"), width = 10, height = 7.5)
  print(km_plot); dev.off()
} else {
  message("package 'survminer' not installed: install.packages('survminer') for the KM plot with risk table; using base graphics")
  pdf(paste0(out_prefix, "_KM_by_creditscore.pdf"), width = 9, height = 6)
  plot(survfit(Surv(time, status) ~ CreditScoreRange, data = loans_data), col = 1:nlevels(loans_data$CreditScoreRange), lwd = 2,
       xlab = "Time (Days)", ylab = "Survival Probability", main = "Kaplan-Meier Survival Curves stratified by Credit Score")
  legend("bottomleft", levels(loans_data$CreditScoreRange), col = 1:nlevels(loans_data$CreditScoreRange), lwd = 2, bty = "n"); dev.off()
}

# ============================================================
# FIT BOTH MODELS ON ALL ANALYSED LOANS
# ============================================================
cat("\nFitting the MIXTURE cure model (design:", design, ")\n")
fit_mix <- fit_mixture(loans_data)
print(fit_mix)
cat("\nFitting the BOX-COX transformation cure model (design:", design, ")\n")
fit_bct <- fit_boxcox(loans_data)
print(fit_bct)

for (nm_ in c("mix", "bct")) {
  ft <- get(paste0("fit_", nm_))
  cat(sprintf("\n[%s] cure fraction quantiles:\n", nm_)); print(quantile(ft$cure_fraction, c(0.05, 0.25, 0.5, 0.75, 0.95)))
  cat(sprintf("[%s] mean cure fraction by credit-score group:\n", nm_)); print(tapply(ft$cure_fraction, loans_data$CreditScoreRange, mean))
}
km_end <- summary(survfit(Surv(time, status) ~ 1, data = loans_data), times = max(loans_data$time))$surv
cat(sprintf("Kaplan-Meier survival at the end of follow-up: %.3f (a cure fraction above this value is not supported by the data)\n", km_end))
if (fit_bct$alpha > 0.99)
  cat("alpha-hat is at the upper boundary: the BCT model coincides with the logistic mixture cure model (cure coefficients equal up to sign).\n")

# estimates side by side (cure coefficients: mixture = log-odds of cure, BCT = eta, opposite signs when alpha = 1)
est_all <- rbind(cbind(Model = "Mixture cure (flexsurvcure)", fit_mix$results),
                 cbind(Model = "Box-Cox cure", fit_bct$results))
write.csv(est_all, paste0(out_prefix, "_estimates.csv"), row.names = FALSE)

# model comparison table including the no-cure benchmarks
cmp <- data.frame(Model = c("Mixture cure (lognormal)", "Box-Cox cure (lognormal)"),
                  LogLik = c(fit_mix$logLik, fit_bct$logLik), npar = c(fit_mix$n_par, fit_bct$n_par),
                  AIC = c(fit_mix$AIC, fit_bct$AIC), BIC = c(fit_mix$BIC, fit_bct$BIC))
if (!is.null(results_std)) { results_std$Model <- paste0("no-cure ", results_std$Model); cmp <- rbind(cmp, results_std[, names(cmp)]) }
cmp <- cmp[order(cmp$AIC), ]
cat("\n========== ALL MODELS, ordered by AIC (same data, maximum likelihood) ==========\n"); print(cmp, digits = 6, row.names = FALSE)
write.csv(cmp, paste0(out_prefix, "_model_comparison.csv"), row.names = FALSE)

# ============================================================
# HOLD-OUT SPLIT (70/30, stratified on status); both models refitted on the training part
# ============================================================
set.seed(2000)
idx_train <- unlist(lapply(split(seq_len(nrow(loans_data)), loans_data$status),
                           function(i) sample(i, size = round(0.7 * length(i)))))
loans_train <- loans_data[sort(idx_train), ]
loans_test  <- loans_data[-idx_train, ]
cat(sprintf("\nTrain n = %d (%d defaults, %.1f%%) | test n = %d (%d defaults, %.1f%%)\n",
            nrow(loans_train), sum(loans_train$status), 100 * mean(loans_train$status),
            nrow(loans_test), sum(loans_test$status), 100 * mean(loans_test$status)))
cat("\n---- training fits ----\n")
fit_mix_train <- fit_mixture(loans_train, verbose = FALSE); print(fit_mix_train)
fit_bct_train <- fit_boxcox(loans_train, verbose = FALSE); print(fit_bct_train)

horizon <- median(loans_train$time[loans_train$status == 1])
for (nm_ in c("mixture", "boxcox")) {
  pr <- if (nm_ == "mixture") predict_mixture_survival(fit_mix_train, loans_test, horizon) else predict_boxcox_survival(fit_bct_train, loans_test, horizon)
  cidx <- survival::concordance(Surv(time, status) ~ I(1 - pr[, 1]), data = loans_test, reverse = TRUE)
  cat(sprintf("Hold-out concordance, %s (risk = 1 - S(%.0f days)): %.3f (SE %.3f)\n", nm_, horizon, cidx$concordance, sqrt(cidx$var)))
}

# calibration of the marginal curve on the test set (PDF)
tg <- seq(0, max(loans_test$time), length.out = 60)
S_mix <- colMeans(predict_mixture_survival(fit_mix_train, loans_test, tg))
S_bct <- colMeans(predict_boxcox_survival(fit_bct_train, loans_test, tg))
km <- survfit(Surv(time, status) ~ 1, data = loans_test)
pdf(paste0(out_prefix, "_calibration.pdf"), width = 9, height = 6)
plot(km, conf.int = FALSE, xlab = "Days", ylab = "Survival (no default)", main = "Test set: Kaplan-Meier vs fitted marginal survival")
lines(tg, S_mix, lty = 2, col = "red", lwd = 2)
lines(tg, S_bct, lty = 3, col = "blue", lwd = 2)
legend("bottomleft", c("Kaplan-Meier", "Mixture cure", "Box-Cox cure"), lty = c(1, 2, 3), col = c("black", "red", "blue"), lwd = 2, bty = "n")
dev.off()

# ============================================================
# SurvSHAP(t) -- MIXTURE CURE MODEL ONLY (survex API: check the installed version if this block fails)
# ============================================================
if (run_shap) {
  library(survex)
  vars <- c("DebtToIncomeRatio", if (nlevels(loans_data$Term) > 1) "Term", "BorrowerRate", "CreditScoreRange",
            "BankcardUtilization", "InquiriesLast6Months", "OpenCreditLines", "IsBorrowerHomeowner")
  X_train <- as.data.frame(loans_train[, vars]); X_test <- as.data.frame(loans_test[, vars])
  time_points <- if (shap_times == "unique") sort(unique(loans_test$time[loans_test$time > 0])) else seq(1, max(loans_test$time), length.out = 100)
  cat("SurvSHAP(t) on", length(time_points), "time points\n")
  predict_risk_fun <- function(model, newdata) 1 - predict_mixture_survival(model, newdata, horizon)[, 1]
  explainer <- survex::explain_survival(
    model = fit_mix_train, data = X_train, y = Surv(loans_train$time, loans_train$status),
    predict_survival_function = predict_mixture_survival, predict_risk_function = predict_risk_fun,
    times = time_points, label = "Mixture cure model (lognormal)")
  set.seed(1)
  shap_obs <- X_test[sample(nrow(X_test), min(shap_n, nrow(X_test))), , drop = FALSE]
  cat("Explaining", nrow(shap_obs), "test loans\n")
  shap_result <- tryCatch(
    survex::model_survshap(explainer, new_observation = shap_obs, N = 100,
                           calculation_method = "kernelshap", aggregation_method = "integral"),
    error = function(e) { message("SurvSHAP failed: ", conditionMessage(e)); NULL })
  if (!is.null(shap_result)) {
    p1 <- plot(shap_result, geom = "importance", max_vars = 12) + ggtitle("Mixture cure model: global feature importance over time") + theme_minimal()
    p2 <- plot(shap_result, geom = "beeswarm", max_vars = 12) + ggtitle("Mixture cure model: SurvSHAP(t) individual effects") + theme_minimal()
    ggsave(paste0(out_prefix, "_mixture_SHAP_importance.pdf"), p1, width = 12, height = 10)
    ggsave(paste0(out_prefix, "_mixture_SHAP_beeswarm.pdf"), p2, width = 12, height = 10)
  }
}
cat("\nDone. Outputs start with:", out_prefix, "\n")