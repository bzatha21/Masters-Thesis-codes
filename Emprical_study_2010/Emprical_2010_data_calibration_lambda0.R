#!/usr/bin/env Rscript

# ============================================================
# Empirical Phase M2: 2010 SPX model comparison, lambda_skew = 0
# ------------------------------------------------------------
# File name:
#   PhaseM_2010_rheston_calibration_lambda0.R
#
# Purpose:
#   Calibrate the validated rough Heston CF/Lewis pricing engine
#   to the cleaned 2010 SPX option panel from PhaseM_2010_market_panel.R,
#   and compare it with two benchmark mechanisms for skew:
#
#     1. classical Heston: stochastic-volatility benchmark, alpha fixed at 1;
#     2. Variance Gamma: pure-jump exponential Levy benchmark.
#
# Model hierarchy:
#   - rough Heston is the main model;
#   - classical Heston and VG are benchmark models.
#
# Interpretation:
#   This is a one-date empirical calibration case study. It reports
#   IV-surface fit, price errors, ATM-IV errors, ATM-skew errors,
#   maturity-level diagnostics, and model limitations. The script does
#   not force rough Heston to dominate every metric.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(2712010)

# ============================================================
# 0. Run controls and paths
# ============================================================

FAST_TEST <- identical(Sys.getenv("PHASEM_FAST_TEST"), "1")
# Extended robustness is the default for non-fast runs. Set
# PHASEM_EXTENDED_ROBUSTNESS=0 only when a lighter diagnostic run is required.
EXTENDED_ROBUSTNESS <- (!FAST_TEST) && !identical(Sys.getenv("PHASEM_EXTENDED_ROBUSTNESS", "1"), "0")
CLEAN_OUTPUT_DIR <- TRUE
PLOT_FORMAT <- "pdf"
SHOW_PLOTS <- FALSE
SAVE_PLOTS <- TRUE

# Numerical pricing controls.
# Defaults are practical terminal-run settings. They are deliberately
# lighter than the full simulation-validation grid, but the same Lewis
# convention and fractional-Riccati implementation are used.
DU_TARGET <- if (FAST_TEST) 0.125 else as.numeric(Sys.getenv("PHASEM_DU_TARGET", "0.125"))
RICCATI_CORRECTOR_ITERS <- 2
# Default weight on the maturity-wise ATM-skew penalty. The environment
# variable remains available for controlled reruns.
LAMBDA_SKEW <- as.numeric(Sys.getenv("PHASEM_LAMBDA_SKEW", "0"))
OPT_TRACE <- identical(Sys.getenv("PHASEM_OPT_TRACE", "1"), "1")

u_max_for_tau_days <- function(td) {
  if (FAST_TEST) {
    if (td <= 14) return(120)
    if (td <= 45) return(100)
    return(80)
  }
  if (td <= 14) return(180)
  if (td <= 45) return(150)
  if (td <= 100) return(130)
  if (td <= 180) return(120)
  return(110)
}

n_time_for_tau_days <- function(td) {
  if (FAST_TEST) {
    if (td <= 14) return(120)
    if (td <= 45) return(100)
    return(80)
  }
  if (td <= 14) return(360)
  if (td <= 45) return(300)
  if (td <= 100) return(260)
  return(220)
}

MAX_CALIB_EXPIRIES <- if (FAST_TEST) 3L else as.integer(Sys.getenv("PHASEM_MAX_CALIB_EXPIRIES", "8"))
MAX_QUOTES_PER_EXPIRY <- if (FAST_TEST) 5L else as.integer(Sys.getenv("PHASEM_MAX_QUOTES_PER_EXPIRY", "15"))
MAXIT_ROUGH <- if (FAST_TEST) 80L else as.integer(Sys.getenv("PHASEM_MAXIT_ROUGH", if (EXTENDED_ROBUSTNESS) "250" else "120"))
MAXIT_HESTON <- if (FAST_TEST) 70L else as.integer(Sys.getenv("PHASEM_MAXIT_HESTON", if (EXTENDED_ROBUSTNESS) "220" else "100"))
MAXIT_VG <- if (FAST_TEST) 90L else as.integer(Sys.getenv("PHASEM_MAXIT_VG", if (EXTENDED_ROBUSTNESS) "250" else "180"))

base_dir <- Sys.getenv("R_THESIS_BASE")
if (identical(base_dir, "")) {
  if (dir.exists("/work/Home/R_Thesis")) {
    base_dir <- "/work/Home/R_Thesis"
  } else if (dir.exists(path.expand("~/R_Thesis"))) {
    base_dir <- path.expand("~/R_Thesis")
  } else {
    base_dir <- getwd()
  }
}

in_dir <- file.path(base_dir, "data", "EmpiricalStudy", "phaseM_2010_market_panel_output")
out_dir <- file.path(base_dir, "data", "EmpiricalStudy", "phaseM_2010_rheston_calibration_lambda0_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_OUTPUT_DIR) {
  old_files <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)
}

cat("Input market-panel directory:\n", in_dir, "\n\n", sep = "")
cat("Output calibration directory:\n", out_dir, "\n\n", sep = "")
cat("FAST_TEST mode:", FAST_TEST, "\n")
cat("EXTENDED_ROBUSTNESS mode:", EXTENDED_ROBUSTNESS, "\n\n")

required_files <- c(
  "phaseM2010_market_calibration_panel.csv",
  "phaseM2010_expiry_summary.csv",
  "phaseM2010_market_atm_skew.csv",
  "phaseM2010_market_atm_skew_loglog_summary.csv",
  "phaseM2010_market_atm_skew_extraction_audit.csv",
  "phaseM2010_market_H_selection_audit.csv",
  "phaseM2010_empirical_H_for_calibration.csv"
)
missing_files <- required_files[!file.exists(file.path(in_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("Missing required Phase M1 files:\n",
       paste(missing_files, collapse = "\n"),
       "\nRun PhaseM_2010_market_panel.R first.")
}

plot_save <- function(plot_obj, stem, width, height, dpi = 180) {
  if (isTRUE(SHOW_PLOTS)) print(plot_obj)
  if (isTRUE(SAVE_PLOTS)) {
    path <- file.path(out_dir, paste0(stem, ".", PLOT_FORMAT))
    if (PLOT_FORMAT == "png") {
      ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
    } else {
      ggsave(path, plot_obj, width = width, height = height)
    }
  }
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  sqrt(mean(x^2))
}

safe_mae <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(abs(x))
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

# ============================================================
# 1. Black--Scholes functions
# ============================================================

black_call_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) {
    return(D * max(F - K, 0))
  }
  vol_sqrt <- sigma * sqrt(tau)
  d1 <- (log(F / K) + 0.5 * sigma^2 * tau) / vol_sqrt
  d2 <- d1 - vol_sqrt
  D * (F * pnorm(d1) - K * pnorm(d2))
}

black_vega_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) return(0)
  vol_sqrt <- sigma * sqrt(tau)
  d1 <- (log(F / K) + 0.5 * sigma^2 * tau) / vol_sqrt
  D * F * dnorm(d1) * sqrt(tau)
}

implied_vol_black_call <- function(price, F, K, tau, D = 1, tol = 1e-8) {
  intrinsic <- D * max(F - K, 0)
  upper <- D * F
  if (!is.finite(price) || !is.finite(F) || !is.finite(K) || !is.finite(tau) || !is.finite(D)) return(NA_real_)
  if (tau <= 0 || F <= 0 || K <= 0 || D <= 0) return(NA_real_)
  if (price < intrinsic - 1e-8 || price > upper + 1e-8) return(NA_real_)
  price <- min(max(price, intrinsic + 1e-12), upper - 1e-12)
  f <- function(sig) black_call_forward(F, K, tau, sig, D) - price
  lo <- 1e-5
  hi <- 5.0
  flo <- f(lo)
  fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  uniroot(f, lower = lo, upper = hi, tol = tol)$root
}

# ============================================================
# 2. Lewis control and rough Heston CF
# ============================================================

black_cf_log_forward <- function(z, tau, sigma) {
  exp(-0.5 * sigma^2 * tau * (z^2 + 1i * z))
}

lewis_call_from_cf_generic <- function(k, tau, cf_vec, u_grid, F, D) {
  K <- F * exp(k)
  du <- u_grid[2] - u_grid[1]
  w <- rep(du, length(u_grid))
  w[1L] <- 0.5 * du
  w[length(w)] <- 0.5 * du
  vals <- exp(-1i * u_grid * k) * cf_vec / (u_grid^2 + 0.25)
  integral <- sum(w * Re(vals))
  price <- D * (F - sqrt(F * K) / pi * integral)
  Re(price)
}

is_finite_complex <- function(z) {
  is.finite(Re(z)) && is.finite(Im(z))
}

make_abm_precomp <- function(T, n_time, alpha) {
  dt <- T / n_time
  n <- n_time
  lag <- 1:n
  
  b_lag <- (dt^alpha / base::gamma(alpha + 1)) *
    (lag^alpha - (lag - 1)^alpha)
  
  a_list <- vector("list", n)
  for (k in 1:n) {
    a <- numeric(k)
    a[1] <- (k - 1)^(alpha + 1) - (k - alpha - 1) * k^alpha
    if (k >= 2) {
      for (j in 1:(k - 1)) {
        lag_kj <- k - j
        a[j + 1] <- (lag_kj + 1)^(alpha + 1) -
          2 * lag_kj^(alpha + 1) +
          (lag_kj - 1)^(alpha + 1)
      }
    }
    a_list[[k]] <- (dt^alpha / base::gamma(alpha + 2)) * a
  }
  
  beta <- 1 - alpha
  if (abs(beta) < 1e-12) {
    weights_beta <- NULL
  } else {
    lag_beta <- n:1
    weights_beta <- ((lag_beta^beta - (lag_beta - 1)^beta) * dt^beta) /
      base::gamma(beta + 1)
  }
  
  list(
    T = T,
    n_time = n_time,
    alpha = alpha,
    dt = dt,
    b_lag = b_lag,
    a_list = a_list,
    weights_beta = weights_beta
  )
}

rough_heston_cf_one_abm <- function(a, pre, par, corrector_iters = 2) {
  n <- pre$n_time
  dt <- pre$dt
  
  gamma_p <- par$gamma
  theta_p <- par$theta
  nu_p <- par$nu
  rho_p <- par$rho
  V0_p <- par$V0
  alpha_p <- par$alpha
  
  ccoef <- gamma_p * nu_p
  
  rhs_fun <- function(y) {
    0.5 * (-a^2 - 1i * a) +
      gamma_p * (1i * a * rho_p * nu_p - 1) * y +
      0.5 * ccoef^2 * y^2
  }
  
  h <- complex(length = n + 1L)
  fvals <- complex(length = n + 1L)
  h[1L] <- 0 + 0i
  fvals[1L] <- rhs_fun(h[1L])
  
  ok <- TRUE
  fail_step <- NA_integer_
  
  for (k in 1:n) {
    h_pred <- sum(pre$b_lag[k:1] * fvals[seq_len(k)])
    if (!is_finite_complex(h_pred) || Mod(h_pred) > 1e10) {
      ok <- FALSE
      fail_step <- k
      break
    }
    
    a_weights <- pre$a_list[[k]]
    base_term <- sum(a_weights * fvals[seq_len(k)])
    h_new <- h_pred
    
    for (it in seq_len(corrector_iters)) {
      f_new <- rhs_fun(h_new)
      h_new <- base_term +
        (dt^pre$alpha / base::gamma(pre$alpha + 2)) * f_new
      if (!is_finite_complex(h_new) || Mod(h_new) > 1e10) {
        ok <- FALSE
        fail_step <- k
        break
      }
    }
    if (!ok) break
    
    h[k + 1L] <- h_new
    fvals[k + 1L] <- rhs_fun(h_new)
    if (!is_finite_complex(fvals[k + 1L]) || Mod(fvals[k + 1L]) > 1e12) {
      ok <- FALSE
      fail_step <- k
      break
    }
  }
  
  if (!ok) {
    return(list(cf = NA_complex_, ok = FALSE, fail_step = fail_step, max_mod_h = safe_max(Mod(h))))
  }
  
  int_h <- dt * (sum(h) - 0.5 * h[1L] - 0.5 * h[length(h)])
  I_beta_h_T <- if (abs(1 - alpha_p) < 1e-12) {
    h[length(h)]
  } else {
    sum(pre$weights_beta * h[seq_len(n)])
  }
  
  exponent <- gamma_p * theta_p * int_h + V0_p * I_beta_h_T
  if (!is_finite_complex(exponent) || Re(exponent) > 700) {
    return(list(cf = NA_complex_, ok = FALSE, fail_step = NA_integer_, max_mod_h = safe_max(Mod(h))))
  }
  
  cf <- exp(exponent)
  list(cf = cf, ok = is_finite_complex(cf), fail_step = NA_integer_, max_mod_h = safe_max(Mod(h)))
}

compute_cf_grid_for_maturity <- function(T, u_grid, par, n_time, corrector_iters = 2) {
  pre <- make_abm_precomp(T = T, n_time = n_time, alpha = par$alpha)
  cf <- complex(length = length(u_grid))
  ok <- logical(length(u_grid))
  max_mod_h <- numeric(length(u_grid))
  
  for (idx in seq_along(u_grid)) {
    a <- u_grid[idx] - 0.5i
    res <- rough_heston_cf_one_abm(a, pre, par, corrector_iters)
    cf[idx] <- res$cf
    ok[idx] <- isTRUE(res$ok)
    max_mod_h[idx] <- res$max_mod_h
  }
  
  tibble(
    u = u_grid,
    cf_real = Re(cf),
    cf_imag = Im(cf),
    cf_ok = ok,
    max_mod_h = max_mod_h
  )
}

price_model_panel <- function(panel, par, corrector_iters = 2) {
  out_list <- list()
  idx <- 1L
  bad_cf_total <- 0L
  clipped_total <- 0L
  
  for (td in sort(unique(panel$tau_days))) {
    d <- panel %>% filter(tau_days == td) %>% arrange(k)
    Tval <- unique(d$tau)[1]
    Fval <- unique(d$F)[1]
    Dval <- unique(d$D)[1]
    
    u_max <- u_max_for_tau_days(td)
    n_u <- ceiling(u_max / DU_TARGET)
    u_grid <- seq(0, u_max, length.out = n_u + 1L)
    n_time <- n_time_for_tau_days(td)
    
    cf_tbl <- compute_cf_grid_for_maturity(
      T = Tval,
      u_grid = u_grid,
      par = par,
      n_time = n_time,
      corrector_iters = corrector_iters
    )
    
    bad_cf <- sum(!cf_tbl$cf_ok | !is.finite(cf_tbl$cf_real) | !is.finite(cf_tbl$cf_imag))
    bad_cf_total <- bad_cf_total + bad_cf
    
    if (bad_cf > 0L) {
      return(list(ok = FALSE, panel = NULL, bad_cf_total = bad_cf_total, clipped_total = clipped_total))
    }
    
    cf_vec <- cf_tbl$cf_real + 1i * cf_tbl$cf_imag
    
    rows <- lapply(seq_len(nrow(d)), function(i) {
      kval <- d$k[i]
      Kval <- d$K[i]
      price_raw <- lewis_call_from_cf_generic(
        k = kval,
        tau = Tval,
        cf_vec = cf_vec,
        u_grid = u_grid,
        F = Fval,
        D = Dval
      )
      intrinsic <- Dval * max(Fval - Kval, 0)
      upper <- Dval * Fval
      price <- min(max(price_raw, intrinsic + 1e-10), upper - 1e-10)
      clipped <- abs(price - price_raw) > 1e-7
      iv <- implied_vol_black_call(price, Fval, Kval, Tval, Dval)
      
      tibble(
        row_id = d$row_id[i],
        model_price_raw = price_raw,
        model_price = price,
        model_iv = iv,
        model_price_clipped = clipped,
        model_n_time = n_time,
        model_u_max = u_max,
        model_n_u = n_u
      )
    })
    
    out_list[[idx]] <- bind_rows(rows)
    clipped_total <- clipped_total + sum(out_list[[idx]]$model_price_clipped, na.rm = TRUE)
    idx <- idx + 1L
  }
  
  fit_panel <- panel %>% left_join(bind_rows(out_list), by = "row_id")
  list(ok = TRUE, panel = fit_panel, bad_cf_total = bad_cf_total, clipped_total = clipped_total)
}


price_rheston_panel <- price_model_panel

# ============================================================
# 2b. Variance Gamma benchmark CF and Lewis pricing
# ============================================================
# Exponential Levy VG benchmark under forward normalization.
# The raw VG increment X_T has characteristic function
#   phi_X(z) = (1 - i theta kappa z + 0.5 sigma^2 kappa z^2)^(-T/kappa).
# The log-forward variable is X_T - T log E[exp(X_1)], so that phi(-i)=1.

make_vg_par <- function(sigma, theta_vg, kappa) {
  list(sigma = sigma, theta_vg = theta_vg, kappa = kappa)
}

vg_valid_par <- function(par) {
  is.finite(par$sigma) && is.finite(par$theta_vg) && is.finite(par$kappa) &&
    par$sigma > 0 && par$kappa > 0 &&
    is.finite(1 - par$theta_vg * par$kappa - 0.5 * par$sigma^2 * par$kappa) &&
    (1 - par$theta_vg * par$kappa - 0.5 * par$sigma^2 * par$kappa) > 1e-10
}

vg_cf_log_forward <- function(z, tau, par) {
  if (!vg_valid_par(par)) {
    return(rep(NA_complex_, length(z)))
  }
  
  sigma <- par$sigma
  theta <- par$theta_vg
  kappa <- par$kappa
  
  den_mart <- 1 - theta * kappa - 0.5 * sigma^2 * kappa
  log_m1 <- -(1 / kappa) * log(den_mart)
  
  base <- 1 - 1i * theta * kappa * z + 0.5 * sigma^2 * kappa * z^2
  cf_raw <- exp((-tau / kappa) * log(base))
  
  # Forward normalization. This enforces E[exp(log-forward)] = 1.
  exp(-1i * z * tau * log_m1) * cf_raw
}

price_vg_panel <- function(panel, par, corrector_iters = NULL) {
  out_list <- list()
  idx <- 1L
  clipped_total <- 0L
  
  if (!vg_valid_par(par)) {
    return(list(ok = FALSE, panel = NULL, bad_cf_total = NA_integer_, clipped_total = NA_integer_))
  }
  
  for (td in sort(unique(panel$tau_days))) {
    d <- panel %>% filter(tau_days == td) %>% arrange(k)
    Tval <- unique(d$tau)[1]
    Fval <- unique(d$F)[1]
    Dval <- unique(d$D)[1]
    
    u_max <- u_max_for_tau_days(td)
    n_u <- ceiling(u_max / DU_TARGET)
    u_grid <- seq(0, u_max, length.out = n_u + 1L)
    
    cf_vec <- vg_cf_log_forward(u_grid - 0.5i, Tval, par)
    bad_cf <- sum(!is.finite(Re(cf_vec)) | !is.finite(Im(cf_vec)))
    if (bad_cf > 0L) {
      return(list(ok = FALSE, panel = NULL, bad_cf_total = bad_cf, clipped_total = clipped_total))
    }
    
    rows <- lapply(seq_len(nrow(d)), function(i) {
      kval <- d$k[i]
      Kval <- d$K[i]
      price_raw <- lewis_call_from_cf_generic(
        k = kval,
        tau = Tval,
        cf_vec = cf_vec,
        u_grid = u_grid,
        F = Fval,
        D = Dval
      )
      intrinsic <- Dval * max(Fval - Kval, 0)
      upper <- Dval * Fval
      price <- min(max(price_raw, intrinsic + 1e-10), upper - 1e-10)
      clipped <- abs(price - price_raw) > 1e-7
      iv <- implied_vol_black_call(price, Fval, Kval, Tval, Dval)
      
      tibble(
        row_id = d$row_id[i],
        model_price_raw = price_raw,
        model_price = price,
        model_iv = iv,
        model_price_clipped = clipped,
        model_n_time = NA_integer_,
        model_u_max = u_max,
        model_n_u = n_u
      )
    })
    
    out_list[[idx]] <- bind_rows(rows)
    clipped_total <- clipped_total + sum(out_list[[idx]]$model_price_clipped, na.rm = TRUE)
    idx <- idx + 1L
  }
  
  fit_panel <- panel %>% left_join(bind_rows(out_list), by = "row_id")
  list(ok = TRUE, panel = fit_panel, bad_cf_total = 0L, clipped_total = clipped_total)
}

# ============================================================
# 3. Read market panel and construct calibration subset
# ============================================================

panel_all <- read_csv(file.path(in_dir, "phaseM2010_market_calibration_panel.csv"), show_col_types = FALSE)
expiry_summary <- read_csv(file.path(in_dir, "phaseM2010_expiry_summary.csv"), show_col_types = FALSE)
market_atm <- read_csv(file.path(in_dir, "phaseM2010_market_atm_skew.csv"), show_col_types = FALSE)
market_H_summary <- read_csv(file.path(in_dir, "phaseM2010_market_atm_skew_loglog_summary.csv"), show_col_types = FALSE)
market_skew_audit <- read_csv(file.path(in_dir, "phaseM2010_market_atm_skew_extraction_audit.csv"), show_col_types = FALSE)
market_H_selection_audit <- read_csv(file.path(in_dir, "phaseM2010_market_H_selection_audit.csv"), show_col_types = FALSE)
empirical_H_input <- read_csv(file.path(in_dir, "phaseM2010_empirical_H_for_calibration.csv"), show_col_types = FALSE)

# Roughness input from Phase M1: H is estimated from the short-end
# ATM-skew term structure under explicit support and robustness rules.
# Phase M2 consumes that audited decision rather than selecting an
# ambiguous row from a generic log-log summary.
H_decision <- empirical_H_input %>%
  filter(selected_for_calibration, is.finite(H_for_calibration), is.finite(alpha_for_calibration)) %>%
  slice(1L) %>%
  transmute(
    H_source = estimator,
    n_maturities,
    maturity_days_used,
    beta_hat,
    H_hat,
    H_for_calibration,
    alpha_for_calibration,
    H_was_clipped,
    r_squared,
    selection_rule,
    calibration_expiries_days,
    interpretation
  )

if (nrow(H_decision) != 1L) {
  stop("No selected empirical H estimate found in phaseM2010_empirical_H_for_calibration.csv")
}

parse_day_list <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(integer(0))
  as.integer(trimws(unlist(strsplit(as.character(x), ","))))
}

HEADLINE_H_DAYS <- parse_day_list(H_decision$maturity_days_used[1])

H_FIXED <- as.numeric(H_decision$H_for_calibration[1])
ALPHA_FIXED <- as.numeric(H_decision$alpha_for_calibration[1])
if (!is.finite(H_FIXED) || !is.finite(ALPHA_FIXED) || ALPHA_FIXED <= 0.5 || ALPHA_FIXED >= 1.0) {
  stop("Invalid fixed rough Heston alpha from Phase M1 H estimate: ", ALPHA_FIXED)
}

write_csv(H_decision, file.path(out_dir, "phaseM2010_rough_H_fixed_from_market_skew.csv"))
cat("Fixed rough Heston H from audited market ATM-skew:", H_FIXED, "\n")
cat("Fixed rough Heston alpha:", ALPHA_FIXED, "\n")
cat("H source:", H_decision$H_source[1], "\n")
cat("H maturities:", H_decision$maturity_days_used[1], "\n")
cat("Headline H-day vector:", paste(HEADLINE_H_DAYS, collapse = ", "), "\n\n")



panel0 <- panel_all %>%
  mutate(
    row_id = row_number(),
    tau_days = as.integer(tau_days),
    tau = as.numeric(tau),
    k = as.numeric(k),
    K = as.numeric(K),
    F = as.numeric(F),
    D = as.numeric(D),
    iv_mkt = as.numeric(iv_mkt),
    weight_iv = as.numeric(weight_iv),
    calibration_flag = as.logical(calibration_flag),
    holdout_flag = as.logical(holdout_flag)
  ) %>%
  filter(
    tau_days <= 365,
    abs(k) <= 0.35,
    is.finite(iv_mkt),
    is.finite(weight_iv),
    weight_iv > 0
  )

selected_expiries <- panel0 %>%
  filter(calibration_flag) %>%
  distinct(tau_days) %>%
  arrange(tau_days) %>%
  slice_head(n = MAX_CALIB_EXPIRIES) %>%
  pull(tau_days)

# IV calibration remains capped for speed. The ATM-skew object is not capped.
SKEW_WINDOW <- as.numeric(Sys.getenv("PHASEM_SKEW_WINDOW", "0.055"))
SKEW_OBJECTIVE_MODE <- Sys.getenv("PHASEM_SKEW_OBJECTIVE_MODE", "headline")

if (!is.finite(SKEW_WINDOW) || SKEW_WINDOW <= 0) {
  stop("Invalid PHASEM_SKEW_WINDOW: ", SKEW_WINDOW)
}

SKEW_OBJECTIVE_DAYS <- if (identical(SKEW_OBJECTIVE_MODE, "all")) {
  selected_expiries
} else if (identical(SKEW_OBJECTIVE_MODE, "headline")) {
  intersect(HEADLINE_H_DAYS, selected_expiries)
} else {
  stop("Unknown PHASEM_SKEW_OBJECTIVE_MODE. Use 'headline' or 'all'.")
}

if (length(SKEW_OBJECTIVE_DAYS) < 3L) {
  stop(
    "Too few maturities in SKEW_OBJECTIVE_DAYS: ",
    paste(SKEW_OBJECTIVE_DAYS, collapse = ", ")
  )
}

panel_work <- panel0 %>%
  filter(tau_days %in% selected_expiries) %>%
  group_by(tau_days) %>%
  arrange(abs(k), .by_group = TRUE) %>%
  mutate(rank_near_atm = row_number()) %>%
  ungroup() %>%
  filter(rank_near_atm <= MAX_QUOTES_PER_EXPIRY | holdout_flag) %>%
  mutate(
    calib_used = calibration_flag & rank_near_atm <= MAX_QUOTES_PER_EXPIRY,
    holdout_used = holdout_flag & rank_near_atm <= MAX_QUOTES_PER_EXPIRY
  )

calib_panel <- panel_work %>% filter(calib_used)
holdout_panel <- panel_work %>% filter(holdout_used)

# Dedicated full near-ATM panel for skew penalty and skew diagnostics.
skew_panel <- panel0 %>%
  filter(
    tau_days %in% selected_expiries,
    abs(k) <= SKEW_WINDOW,
    is.finite(iv_mkt),
    iv_mkt > 0
  ) %>%
  group_by(tau_days) %>%
  arrange(k, .by_group = TRUE) %>%
  mutate(
    skew_diag_used = TRUE,
    skew_objective_used = tau_days %in% SKEW_OBJECTIVE_DAYS
  ) %>%
  ungroup()

skew_support_summary <- skew_panel %>%
  group_by(tau_days, tau) %>%
  summarise(
    n_skew_local = n(),
    n_skew_left = sum(k < -1e-10),
    n_skew_right = sum(k > 1e-10),
    two_sided_skew_support = n_skew_left >= 1L & n_skew_right >= 1L,
    used_in_skew_objective = any(skew_objective_used),
    .groups = "drop"
  ) %>%
  arrange(tau_days)

bad_skew_support <- skew_support_summary %>%
  filter(used_in_skew_objective, !two_sided_skew_support)

if (nrow(bad_skew_support) > 0L) {
  stop(
    "Skew objective has maturities without two-sided near-ATM support: ",
    paste(bad_skew_support$tau_days, collapse = ", ")
  )
}

calib_ids <- calib_panel$row_id
holdout_ids <- holdout_panel$row_id
skew_diag_ids <- skew_panel$row_id
skew_objective_ids <- skew_panel %>%
  filter(skew_objective_used) %>%
  pull(row_id)

make_pricing_panel <- function(ids) {
  ids <- unique(ids)
  panel0 %>%
    filter(row_id %in% ids) %>%
    mutate(
      calib_used = row_id %in% calib_ids,
      holdout_used = row_id %in% holdout_ids,
      skew_diag_used = row_id %in% skew_diag_ids,
      skew_objective_used = row_id %in% skew_objective_ids
    )
}

objective_panel <- make_pricing_panel(c(calib_ids, skew_objective_ids))
eval_panel <- make_pricing_panel(c(calib_ids, holdout_ids, skew_diag_ids))

write_csv(
  skew_support_summary,
  file.path(out_dir, "phaseM2010_full_near_atm_skew_panel_support.csv")
)

if (nrow(calib_panel) < 10L) {
  stop("Too few calibration quotes after filtering. Inspect PhaseM_2010_market_panel.R outputs.")
}

panel_summary <- tibble(
  fast_test = FAST_TEST,
  extended_robustness = EXTENDED_ROBUSTNESS,
  lambda_skew = LAMBDA_SKEW,
  skew_window = SKEW_WINDOW,
  skew_objective_mode = SKEW_OBJECTIVE_MODE,
  skew_objective_maturities_days = paste(SKEW_OBJECTIVE_DAYS, collapse = ", "),
  n_skew_diagnostic_quotes = length(skew_diag_ids),
  n_skew_objective_quotes = length(skew_objective_ids),
  selected_expiries = paste(selected_expiries, collapse = ", "),
  n_calibration_quotes = nrow(calib_panel),
  n_holdout_quotes = nrow(holdout_panel),
  max_calib_expiries = MAX_CALIB_EXPIRIES,
  max_quotes_per_expiry = MAX_QUOTES_PER_EXPIRY,
  du_target = DU_TARGET,
  fixed_H_from_market_skew = H_FIXED,
  fixed_alpha_from_market_skew = ALPHA_FIXED,
  headline_H_maturities_days = paste(HEADLINE_H_DAYS, collapse = ", "),
  maxit_rough = MAXIT_ROUGH,
  maxit_heston = MAXIT_HESTON,
  maxit_vg = MAXIT_VG
)

write_csv(panel_summary, file.path(out_dir, "phaseM2010_calibration_panel_summary.csv"))
print(panel_summary)

# ============================================================
# 4. Calibration parameter transforms
# ============================================================

logit <- function(x, lo, hi) log((x - lo) / (hi - x))
inv_logit <- function(y, lo, hi) lo + (hi - lo) / (1 + exp(-y))

bounds <- list(
  gamma = c(0.05, 5.00),
  theta = c(0.005, 0.50),
  nu = c(0.03, 1.50),
  rho = c(-0.95, -0.02),
  V0 = c(0.005, 0.50),
  H = c(0.03, 0.49)
)

make_par <- function(gamma, theta, nu, rho, V0, H = 0.10, alpha_fixed = NULL) {
  alpha <- if (is.null(alpha_fixed)) H + 0.5 else alpha_fixed
  list(gamma = gamma, theta = theta, nu = nu, rho = rho, V0 = V0, H = alpha - 0.5, alpha = alpha)
}

pack_rough <- function(p) {
  c(
    logit(p$gamma, bounds$gamma[1], bounds$gamma[2]),
    logit(p$theta, bounds$theta[1], bounds$theta[2]),
    logit(p$nu, bounds$nu[1], bounds$nu[2]),
    logit(p$rho, bounds$rho[1], bounds$rho[2]),
    logit(p$V0, bounds$V0[1], bounds$V0[2])
  )
}

unpack_rough <- function(x) {
  make_par(
    gamma = inv_logit(x[1], bounds$gamma[1], bounds$gamma[2]),
    theta = inv_logit(x[2], bounds$theta[1], bounds$theta[2]),
    nu = inv_logit(x[3], bounds$nu[1], bounds$nu[2]),
    rho = inv_logit(x[4], bounds$rho[1], bounds$rho[2]),
    V0 = inv_logit(x[5], bounds$V0[1], bounds$V0[2]),
    alpha_fixed = ALPHA_FIXED
  )
}

pack_heston <- function(p) {
  c(
    logit(p$gamma, bounds$gamma[1], bounds$gamma[2]),
    logit(p$theta, bounds$theta[1], bounds$theta[2]),
    logit(p$nu, bounds$nu[1], bounds$nu[2]),
    logit(p$rho, bounds$rho[1], bounds$rho[2]),
    logit(p$V0, bounds$V0[1], bounds$V0[2])
  )
}

unpack_heston <- function(x) {
  make_par(
    gamma = inv_logit(x[1], bounds$gamma[1], bounds$gamma[2]),
    theta = inv_logit(x[2], bounds$theta[1], bounds$theta[2]),
    nu = inv_logit(x[3], bounds$nu[1], bounds$nu[2]),
    rho = inv_logit(x[4], bounds$rho[1], bounds$rho[2]),
    V0 = inv_logit(x[5], bounds$V0[1], bounds$V0[2]),
    alpha_fixed = 1.0
  )
}

# Data-informed starting values.
first_atm_iv <- market_atm %>%
  filter(tau_days %in% selected_expiries, is.finite(atm_iv)) %>%
  arrange(tau_days) %>%
  pull(atm_iv)

V0_start <- if (length(first_atm_iv) > 0L) pmin(pmax(first_atm_iv[1]^2, 0.01), 0.20) else 0.045
theta_start <- if (length(first_atm_iv) > 0L) pmin(pmax(tail(first_atm_iv, 1)^2, 0.01), 0.20) else 0.04
atm_sigma_start <- if (length(first_atm_iv) > 0L) pmin(pmax(median(first_atm_iv, na.rm = TRUE), 0.08), 0.50) else 0.20

rough_starts <- list(
  make_par(1.5, 0.040, 0.45, -0.70, 0.045, alpha_fixed = ALPHA_FIXED),
  make_par(0.8, theta_start, 0.35, -0.65, V0_start, alpha_fixed = ALPHA_FIXED),
  make_par(2.2, 0.035, 0.60, -0.80, 0.055, alpha_fixed = ALPHA_FIXED),
  make_par(0.25, 0.3156, 0.331, -0.681, 0.0392, alpha_fixed = ALPHA_FIXED),
  make_par(1.2, 0.060, 0.30, -0.55, V0_start, alpha_fixed = ALPHA_FIXED)
)

heston_starts <- lapply(rough_starts[1:4], function(p) {
  make_par(p$gamma, p$theta, p$nu, p$rho, p$V0, alpha_fixed = 1.0)
})

# VG parameters: sigma, theta_vg, kappa. Negative theta_vg produces left skew.
vg_bounds <- list(
  sigma = c(0.02, 0.90),
  theta_vg = c(-1.50, 0.50),
  kappa = c(0.02, 2.00)
)

pack_vg <- function(p) {
  c(
    logit(p$sigma, vg_bounds$sigma[1], vg_bounds$sigma[2]),
    logit(p$theta_vg, vg_bounds$theta_vg[1], vg_bounds$theta_vg[2]),
    logit(p$kappa, vg_bounds$kappa[1], vg_bounds$kappa[2])
  )
}

unpack_vg <- function(x) {
  make_vg_par(
    sigma = inv_logit(x[1], vg_bounds$sigma[1], vg_bounds$sigma[2]),
    theta_vg = inv_logit(x[2], vg_bounds$theta_vg[1], vg_bounds$theta_vg[2]),
    kappa = inv_logit(x[3], vg_bounds$kappa[1], vg_bounds$kappa[2])
  )
}

vg_starts <- list(
  make_vg_par(sigma = atm_sigma_start, theta_vg = -0.15, kappa = 0.20),
  make_vg_par(sigma = pmax(0.10, 0.85 * atm_sigma_start), theta_vg = -0.30, kappa = 0.45),
  make_vg_par(sigma = pmax(0.12, 1.10 * atm_sigma_start), theta_vg = -0.08, kappa = 0.10),
  make_vg_par(sigma = pmax(0.10, 0.75 * atm_sigma_start), theta_vg = -0.45, kappa = 0.75)
)

if (FAST_TEST) {
  rough_starts <- rough_starts[1:2]
  heston_starts <- heston_starts[1:2]
  vg_starts <- vg_starts[1:2]
} else if (EXTENDED_ROBUSTNESS) {
  # Use all predefined starts for the full multi-start run. This does not
  # certify a global optimum, but it materially reduces dependence on any
  # single local initialisation.
  rough_starts <- rough_starts
  heston_starts <- heston_starts
  vg_starts <- vg_starts
} else {
  rough_starts <- rough_starts[c(1, 2, 4)]
  heston_starts <- heston_starts[c(1, 2, 4)]
  vg_starts <- vg_starts[c(1, 2, 4)]
}

start_design_summary <- tibble(
  fast_test = FAST_TEST,
  extended_robustness = EXTENDED_ROBUSTNESS,
  lambda_skew = LAMBDA_SKEW,
  model = c("rough_heston_fixed_alpha", "classical_heston_alpha_1", "variance_gamma"),
  n_starts = c(length(rough_starts), length(heston_starts), length(vg_starts)),
  maxit = c(MAXIT_ROUGH, MAXIT_HESTON, MAXIT_VG),
  fixed_H_from_market_skew = H_FIXED,
  fixed_alpha_from_market_skew = ALPHA_FIXED
)
write_csv(start_design_summary, file.path(out_dir, "phaseM2010_start_design_summary.csv"))
print(start_design_summary)

# ============================================================
# 5. Numerical control test, objective and calibration
# ============================================================

run_black_lewis_sanity <- function(panel) {
  maturity_rows <- panel %>%
    distinct(tau_days, tau, F, D) %>%
    arrange(tau_days)
  
  test_rows <- tidyr::expand_grid(
    maturity_rows,
    k_test = c(-0.10, 0.00, 0.10)
  )
  
  sigma_test <- 0.22
  
  details <- lapply(seq_len(nrow(test_rows)), function(i) {
    td <- test_rows$tau_days[i]
    Tval <- test_rows$tau[i]
    Fval <- test_rows$F[i]
    Dval <- test_rows$D[i]
    kval <- test_rows$k_test[i]
    Kval <- Fval * exp(kval)
    
    u_max <- u_max_for_tau_days(td)
    n_u <- ceiling(u_max / DU_TARGET)
    u_grid <- seq(0, u_max, length.out = n_u + 1L)
    cf_vec <- sapply(u_grid, function(u) black_cf_log_forward(u - 0.5i, Tval, sigma_test))
    
    p_black <- black_call_forward(Fval, Kval, Tval, sigma_test, Dval)
    p_lewis <- lewis_call_from_cf_generic(kval, Tval, cf_vec, u_grid, Fval, Dval)
    
    tibble(
      tau_days = td,
      tau = Tval,
      k = kval,
      u_max = u_max,
      n_u = n_u,
      price_black = p_black,
      price_lewis = p_lewis,
      abs_error = abs(p_black - p_lewis)
    )
  }) %>% bind_rows()
  
  summary <- details %>%
    summarise(
      max_abs_error = max(abs_error, na.rm = TRUE),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      pass = max_abs_error < 1e-4
    )
  
  list(details = details, summary = summary)
}

black_sanity <- run_black_lewis_sanity(calib_panel)
write_csv(black_sanity$details, file.path(out_dir, "phaseM2010_black_lewis_sanity_details.csv"))
write_csv(black_sanity$summary, file.path(out_dir, "phaseM2010_black_lewis_sanity.csv"))
print(black_sanity$summary)
if (!isTRUE(black_sanity$summary$pass[1])) {
  stop("Black--Scholes/Lewis control failed. Do not calibrate models until the pricing convention is fixed.")
}


local_slope_for_objective <- function(x, y, window = SKEW_WINDOW) {
  ok <- is.finite(x) & is.finite(y) & abs(x) <= window

  d <- tibble(k = x[ok], y = y[ok]) %>%
    arrange(k) %>%
    mutate(local_weight = pmax(0, 1 - abs(k) / window))

  n_left <- sum(d$k < -1e-10)
  n_right <- sum(d$k > 1e-10)

  if (nrow(d) < 3L || n_left < 1L || n_right < 1L) {
    return(NA_real_)
  }

  fit_type <- ifelse(
    nrow(d) >= 5L && n_left >= 2L && n_right >= 2L,
    "local_quadratic",
    "local_linear"
  )

  form <- if (fit_type == "local_quadratic") {
    y ~ k + I(k^2)
  } else {
    y ~ k
  }

  fit <- tryCatch(lm(form, data = d, weights = local_weight), error = function(e) NULL)
  if (is.null(fit) || !"k" %in% names(coef(fit))) {
    return(NA_real_)
  }

  unname(coef(fit)["k"])
}

objective_from_fit <- function(fit_panel) {
  d_iv <- fit_panel %>%
    filter(calib_used) %>%
    mutate(resid = model_iv - iv_mkt) %>%
    filter(is.finite(resid), is.finite(weight_iv))

  if (nrow(d_iv) < 0.8 * nrow(calib_panel)) {
    return(1e7)
  }

  iv_loss <- mean(d_iv$weight_iv * d_iv$resid^2)

  if (!is.finite(LAMBDA_SKEW) || LAMBDA_SKEW <= 0) {
    return(iv_loss)
  }

  d_skew <- fit_panel %>%
    filter(skew_objective_used)

  skew_tbl <- d_skew %>%
    group_by(tau_days) %>%
    summarise(
      market_skew = local_slope_for_objective(k, iv_mkt, window = SKEW_WINDOW),
      model_skew = local_slope_for_objective(k, model_iv, window = SKEW_WINDOW),
      n_skew_local = n(),
      n_left = sum(k < -1e-10),
      n_right = sum(k > 1e-10),
      .groups = "drop"
    ) %>%
    filter(is.finite(market_skew), is.finite(model_skew)) %>%
    mutate(skew_error = model_skew - market_skew)

  if (nrow(skew_tbl) < 3L) {
    return(1e7)
  }

  skew_loss <- mean(skew_tbl$skew_error^2)
  iv_loss + LAMBDA_SKEW * skew_loss
}

calibrate_model <- function(model_name, starts, unpack_fun, pack_fun, price_fun, maxit) {
  cat("\nCalibrating", model_name, "with", length(starts), "starts...\n")
  results <- vector("list", length(starts))
  
  for (sidx in seq_along(starts)) {
    cat("  start", sidx, "of", length(starts), "\n")
    x0 <- pack_fun(starts[[sidx]])
    
    obj <- function(x) {
      par <- unpack_fun(x)
      priced <- price_fun(objective_panel, par, corrector_iters = RICCATI_CORRECTOR_ITERS)
      if (!isTRUE(priced$ok)) return(1e8)
      if (priced$clipped_total > 0L) {
        return(1e7 + 1e5 * priced$clipped_total)
      }
      val <- objective_from_fit(priced$panel)
      if (!is.finite(val)) return(1e8)
      val
    }
    
    fit <- optim(
      par = x0,
      fn = obj,
      method = "Nelder-Mead",
      control = list(
        maxit = maxit,
        reltol = 1e-5,
        trace = ifelse(OPT_TRACE, 1, 0),
        REPORT = 20
      )
    )
    
    par_hat <- unpack_fun(fit$par)
    priced_eval <- price_fun(eval_panel, par_hat, corrector_iters = RICCATI_CORRECTOR_ITERS)
    
    results[[sidx]] <- list(
      fit = fit,
      par = par_hat,
      priced_eval = priced_eval,
      ok = isTRUE(priced_eval$ok)
    )
  }
  
  obj_values <- sapply(results, function(x) x$fit$value)
  best_idx <- which.min(obj_values)
  list(
    model_name = model_name,
    results = results,
    best_idx = best_idx,
    best = results[[best_idx]]
  )
}

rough_cal <- calibrate_model(
  model_name = "rough_heston_fixed_alpha",
  starts = rough_starts,
  unpack_fun = unpack_rough,
  pack_fun = pack_rough,
  price_fun = price_rheston_panel,
  maxit = MAXIT_ROUGH
)

heston_cal <- calibrate_model(
  model_name = "classical_heston_alpha_1",
  starts = heston_starts,
  unpack_fun = unpack_heston,
  pack_fun = pack_heston,
  price_fun = price_rheston_panel,
  maxit = MAXIT_HESTON
)

vg_cal <- calibrate_model(
  model_name = "variance_gamma",
  starts = vg_starts,
  unpack_fun = unpack_vg,
  pack_fun = pack_vg,
  price_fun = price_vg_panel,
  maxit = MAXIT_VG
)

extract_param_tbl <- function(cal) {
  rows <- lapply(seq_along(cal$results), function(i) {
    p <- cal$results[[i]]$par
    
    if (cal$model_name == "variance_gamma") {
      row <- tibble(
        model = cal$model_name,
        start_id = i,
        selected_best = i == cal$best_idx,
        objective = cal$results[[i]]$fit$value,
        convergence = cal$results[[i]]$fit$convergence,
        gamma = NA_real_,
        theta = NA_real_,
        nu = NA_real_,
        rho = NA_real_,
        V0 = NA_real_,
        alpha = NA_real_,
        H = NA_real_,
        vg_sigma = p$sigma,
        vg_theta = p$theta_vg,
        vg_kappa = p$kappa,
        ok_eval = cal$results[[i]]$ok
      )
    } else {
      row <- tibble(
        model = cal$model_name,
        start_id = i,
        selected_best = i == cal$best_idx,
        objective = cal$results[[i]]$fit$value,
        convergence = cal$results[[i]]$fit$convergence,
        gamma = p$gamma,
        theta = p$theta,
        nu = p$nu,
        rho = p$rho,
        V0 = p$V0,
        alpha = p$alpha,
        H = p$H,
        vg_sigma = NA_real_,
        vg_theta = NA_real_,
        vg_kappa = NA_real_,
        ok_eval = cal$results[[i]]$ok
      )
    }
    
    row
  })
  bind_rows(rows)
}

param_tbl <- bind_rows(
  extract_param_tbl(rough_cal),
  extract_param_tbl(heston_cal),
  extract_param_tbl(vg_cal)
)

write_csv(param_tbl, file.path(out_dir, "phaseM2010_calibrated_parameters_all_starts.csv"))
print(param_tbl)

objective_gap_summary <- param_tbl %>%
  group_by(model) %>%
  group_modify(~ {
    d <- .x %>% filter(is.finite(objective)) %>% arrange(objective)
    if (nrow(d) == 0L) {
      return(tibble(
        n_starts = nrow(.x), n_finite_objectives = 0L, n_ok_eval = sum(.x$ok_eval, na.rm = TRUE),
        n_convergence_code_0 = sum(.x$convergence == 0, na.rm = TRUE),
        selected_start_id = NA_integer_, best_objective = NA_real_, second_best_objective = NA_real_,
        objective_gap_second_minus_best = NA_real_, objective_ratio_second_over_best = NA_real_,
        best_convergence_code = NA_integer_
      ))
    }
    best <- d[1, , drop = FALSE]
    second <- if (nrow(d) >= 2L) d[2, , drop = FALSE] else NULL
    second_obj <- if (is.null(second)) NA_real_ else second$objective[1]
    tibble(
      n_starts = nrow(.x),
      n_finite_objectives = nrow(d),
      n_ok_eval = sum(.x$ok_eval, na.rm = TRUE),
      n_convergence_code_0 = sum(.x$convergence == 0, na.rm = TRUE),
      selected_start_id = best$start_id[1],
      best_objective = best$objective[1],
      second_best_objective = second_obj,
      objective_gap_second_minus_best = second_obj - best$objective[1],
      objective_ratio_second_over_best = second_obj / best$objective[1],
      best_convergence_code = best$convergence[1]
    )
  }) %>%
  ungroup() %>%
  mutate(
    robustness_comment = case_when(
      !is.finite(second_best_objective) ~ "Only one finite objective; optimisation robustness is weak.",
      best_convergence_code != 0 ~ "Best start did not report optim convergence code 0; inspect before strong claims.",
      objective_ratio_second_over_best <= 1.15 ~ "Several starts are close to the best objective; local-start sensitivity appears limited.",
      TRUE ~ "Best objective is isolated from the second-best start; treat calibration as local and run/interpret robustness carefully."
    )
  )

write_csv(objective_gap_summary, file.path(out_dir, "phaseM2010_optimization_robustness_summary.csv"))
print(objective_gap_summary)

best_params <- param_tbl %>%
  filter(selected_best) %>%
  arrange(model)

write_csv(best_params, file.path(out_dir, "phaseM2010_calibrated_parameters_best.csv"))

# ============================================================
# 6. Model panels and summaries
# ============================================================

rough_panel <- rough_cal$best$priced_eval$panel %>%
  mutate(model = "rough Heston", model_role = "main", model_family = "rough volatility")

heston_panel <- heston_cal$best$priced_eval$panel %>%
  mutate(model = "classical Heston (alpha=1)", model_role = "benchmark", model_family = "stochastic volatility")

vg_panel <- vg_cal$best$priced_eval$panel %>%
  mutate(model = "Variance Gamma", model_role = "benchmark", model_family = "pure jump Levy")

fit_panel <- bind_rows(rough_panel, heston_panel, vg_panel) %>%
  mutate(
    iv_error = model_iv - iv_mkt,
    price_error = model_price - call_equiv_price,
    sample = case_when(
      calib_used ~ "calibration",
      holdout_used ~ "holdout",
      TRUE ~ "other"
    ),
    moneyness_region = case_when(
      abs(k) <= 0.05 ~ "core_abs_k_le_0p05",
      abs(k) <= 0.10 ~ "near_abs_k_le_0p10",
      TRUE ~ "outer_abs_k_gt_0p10"
    )
  )

error_summary <- fit_panel %>%
  filter(sample %in% c("calibration", "holdout")) %>%
  group_by(model, model_role, model_family, sample) %>%
  summarise(
    n_quotes = n(),
    n_price_clipped = sum(model_price_clipped, na.rm = TRUE),
    n_na_iv = sum(!is.finite(model_iv)),
    iv_rmse = safe_rmse(iv_error),
    iv_mae = safe_mae(iv_error),
    price_rmse = safe_rmse(price_error),
    max_abs_iv_error = safe_max(abs(iv_error)),
    .groups = "drop"
  ) %>%
  arrange(sample, iv_rmse)

error_by_expiry <- fit_panel %>%
  filter(sample %in% c("calibration", "holdout")) %>%
  group_by(model, model_role, model_family, sample, tau_days, tau) %>%
  summarise(
    n_quotes = n(),
    iv_rmse = safe_rmse(iv_error),
    iv_mae = safe_mae(iv_error),
    max_abs_iv_error = safe_max(abs(iv_error)),
    .groups = "drop"
  ) %>%
  arrange(model, sample, tau_days)

error_by_region <- fit_panel %>%
  filter(sample %in% c("calibration", "holdout")) %>%
  group_by(model, model_role, model_family, sample, moneyness_region) %>%
  summarise(
    n_quotes = n(),
    iv_rmse = safe_rmse(iv_error),
    iv_mae = safe_mae(iv_error),
    price_rmse = safe_rmse(price_error),
    max_abs_iv_error = safe_max(abs(iv_error)),
    .groups = "drop"
  ) %>%
  arrange(sample, moneyness_region, iv_rmse)

write_csv(fit_panel, file.path(out_dir, "phaseM2010_model_fit_panel.csv"))
write_csv(error_summary, file.path(out_dir, "phaseM2010_model_error_summary.csv"))
write_csv(error_by_expiry, file.path(out_dir, "phaseM2010_model_error_by_expiry.csv"))
write_csv(error_by_region, file.path(out_dir, "phaseM2010_model_error_by_region.csv"))

print(error_summary)

# ============================================================
# 7. ATM skew comparison
# ============================================================

estimate_atm_skew <- function(df, iv_col, window = 0.055) {
  d <- df %>%
    filter(abs(k) <= window, is.finite(.data[[iv_col]])) %>%
    arrange(k)
  
  n_left <- sum(d$k < -1e-10)
  n_right <- sum(d$k > 1e-10)
  two_sided_support <- n_left >= 1L && n_right >= 1L
  
  if (nrow(d) < 3L || !two_sided_support) {
    return(tibble(
      atm_iv = NA_real_, atm_skew = NA_real_, n_local = nrow(d),
      n_left = n_left, n_right = n_right, two_sided_support = two_sided_support,
      fit_type = "insufficient"
    ))
  }
  
  d <- d %>% mutate(local_weight = pmax(0, 1 - abs(k) / window))
  fit_type <- ifelse(nrow(d) >= 5L && n_left >= 2L && n_right >= 2L, "local_quadratic", "local_linear")
  form <- if (fit_type == "local_quadratic") {
    as.formula(paste0(iv_col, " ~ k + I(k^2)"))
  } else {
    as.formula(paste0(iv_col, " ~ k"))
  }
  
  fit <- tryCatch(lm(form, data = d, weights = local_weight), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      atm_iv = NA_real_, atm_skew = NA_real_, n_local = nrow(d),
      n_left = n_left, n_right = n_right, two_sided_support = two_sided_support,
      fit_type = "failed"
    ))
  }
  
  cf <- coef(fit)
  tibble(
    atm_iv = unname(cf["(Intercept)"]),
    atm_skew = unname(cf["k"]),
    n_local = nrow(d),
    n_left = n_left,
    n_right = n_right,
    two_sided_support = two_sided_support,
    fit_type = fit_type
  )
}

market_skew_eval <- eval_panel %>% filter(skew_diag_used) %>% group_by(tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "iv_mkt", window = 0.055)) %>%
  ungroup() %>%
  mutate(model = "market")

model_skew_eval <- fit_panel %>% filter(skew_diag_used) %>% group_by(model, tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "model_iv", window = 0.055)) %>%
  ungroup()

atm_skew_comparison <- model_skew_eval %>%
  left_join(
    market_skew_eval %>%
      select(tau_days, tau, market_atm_iv = atm_iv, market_atm_skew = atm_skew),
    by = c("tau_days", "tau")
  ) %>%
  mutate(
    atm_iv_error = atm_iv - market_atm_iv,
    atm_skew_error = atm_skew - market_atm_skew,
    abs_atm_skew_error = abs(atm_skew_error),
    in_headline_H_window = tau_days %in% HEADLINE_H_DAYS
  ) %>%
  arrange(model, tau_days)

atm_skew_summary <- atm_skew_comparison %>%
  group_by(model) %>%
  summarise(
    n_expiries = sum(is.finite(atm_skew_error)),
    atm_iv_rmse = safe_rmse(atm_iv_error),
    atm_skew_rmse = safe_rmse(atm_skew_error),
    atm_skew_mae = safe_mae(atm_skew_error),
    max_abs_atm_skew_error = safe_max(abs_atm_skew_error),
    .groups = "drop"
  ) %>%
  arrange(atm_skew_rmse)

atm_skew_summary_headline_H <- atm_skew_comparison %>%
  filter(in_headline_H_window) %>%
  group_by(model) %>%
  summarise(
    n_expiries = sum(is.finite(atm_skew_error)),
    headline_H_maturities_days = paste(sort(unique(tau_days[is.finite(atm_skew_error)])), collapse = ", "),
    atm_iv_rmse = safe_rmse(atm_iv_error),
    atm_skew_rmse = safe_rmse(atm_skew_error),
    atm_skew_mae = safe_mae(atm_skew_error),
    max_abs_atm_skew_error = safe_max(abs_atm_skew_error),
    .groups = "drop"
  ) %>%
  arrange(atm_skew_rmse)

atm_skew_support_audit <- bind_rows(
  market_skew_eval %>%
    mutate(series = "market", model = "market") %>%
    select(series, model, tau_days, tau, atm_iv, atm_skew, n_local, n_left, n_right, two_sided_support, fit_type),
  model_skew_eval %>%
    mutate(series = "model") %>%
    select(series, model, tau_days, tau, atm_iv, atm_skew, n_local, n_left, n_right, two_sided_support, fit_type)
) %>%
  mutate(in_headline_H_window = tau_days %in% HEADLINE_H_DAYS) %>%
  arrange(tau_days, series, model)

write_csv(atm_skew_comparison, file.path(out_dir, "phaseM2010_atm_skew_market_vs_models.csv"))
write_csv(atm_skew_summary, file.path(out_dir, "phaseM2010_atm_skew_error_summary.csv"))
write_csv(atm_skew_summary_headline_H, file.path(out_dir, "phaseM2010_atm_skew_error_summary_headline_H_maturities.csv"))
write_csv(atm_skew_support_audit, file.path(out_dir, "phaseM2010_atm_skew_extraction_support_audit.csv"))

print(atm_skew_summary)
cat("\nATM-skew summary on headline H-selection maturities only:\n")
print(atm_skew_summary_headline_H)

# ============================================================
# 8. Acceptance / interpretation summary
# ============================================================

rheston_martingale_checks <- bind_rows(lapply(sort(unique(eval_panel$tau_days)), function(td) {
  Tval <- unique(eval_panel$tau[eval_panel$tau_days == td])[1]
  nt <- n_time_for_tau_days(td)
  
  check_one <- function(par, model_name) {
    pre <- make_abm_precomp(T = Tval, n_time = nt, alpha = par$alpha)
    cf0 <- rough_heston_cf_one_abm(0 + 0i, pre, par, RICCATI_CORRECTOR_ITERS)
    cf_m1 <- rough_heston_cf_one_abm(-1i, pre, par, RICCATI_CORRECTOR_ITERS)
    tibble(
      model = model_name,
      tau_days = td,
      tau = Tval,
      n_time = nt,
      abs_phi0_minus_1 = Mod(cf0$cf - 1),
      abs_phim1_minus_1 = Mod(cf_m1$cf - 1),
      phi0_ok = isTRUE(cf0$ok),
      phim1_ok = isTRUE(cf_m1$ok)
    )
  }
  
  bind_rows(
    check_one(rough_cal$best$par, "rough Heston"),
    check_one(heston_cal$best$par, "classical Heston (alpha=1)")
  )
}))

vg_martingale_checks <- bind_rows(lapply(sort(unique(eval_panel$tau_days)), function(td) {
  Tval <- unique(eval_panel$tau[eval_panel$tau_days == td])[1]
  cf0 <- vg_cf_log_forward(0 + 0i, Tval, vg_cal$best$par)
  cf_m1 <- vg_cf_log_forward(-1i, Tval, vg_cal$best$par)
  tibble(
    model = "Variance Gamma",
    tau_days = td,
    tau = Tval,
    abs_phi0_minus_1 = Mod(cf0 - 1),
    abs_phim1_minus_1 = Mod(cf_m1 - 1),
    phi0_ok = is.finite(Re(cf0)) && is.finite(Im(cf0)),
    phim1_ok = is.finite(Re(cf_m1)) && is.finite(Im(cf_m1))
  )
}))

martingale_checks <- bind_rows(
  rheston_martingale_checks,
  vg_martingale_checks %>% mutate(n_time = NA_integer_, .before = abs_phi0_minus_1)
)

write_csv(martingale_checks, file.path(out_dir, "phaseM2010_model_martingale_checks.csv"))

ranking_summary <- bind_rows(
  error_summary %>%
    filter(sample == "calibration") %>%
    select(model, model_role, calibration_iv_rmse = iv_rmse),
  error_summary %>%
    filter(sample == "holdout") %>%
    select(model, model_role, holdout_iv_rmse = iv_rmse)
) %>%
  group_by(model, model_role) %>%
  summarise(
    calibration_iv_rmse = safe_max(calibration_iv_rmse),
    holdout_iv_rmse = safe_max(holdout_iv_rmse),
    .groups = "drop"
  ) %>%
  left_join(
    atm_skew_summary %>% select(model, atm_iv_rmse, atm_skew_rmse, max_abs_atm_skew_error),
    by = "model"
  ) %>%
  left_join(
    atm_skew_summary_headline_H %>%
      select(model,
             headline_H_atm_iv_rmse = atm_iv_rmse,
             headline_H_atm_skew_rmse = atm_skew_rmse,
             headline_H_max_abs_atm_skew_error = max_abs_atm_skew_error),
    by = "model"
  ) %>%
  arrange(headline_H_atm_skew_rmse, atm_skew_rmse)

write_csv(ranking_summary, file.path(out_dir, "phaseM2010_model_ranking_summary.csv"))

acceptance <- tibble(
  check = c(
    "Black--Scholes/Lewis control passed",
    "rough Heston calibration completed",
    "classical Heston calibration completed",
    "Variance Gamma calibration completed",
    "rough Heston clipping count on eval panel = 0",
    "classical Heston clipping count on eval panel = 0",
    "Variance Gamma clipping count on eval panel = 0",
    "rough/classical Heston max |phi(0)-1| < 1e-6",
    "rough/classical Heston max |phi(-i)-1| < 1e-6",
    "Variance Gamma max |phi(0)-1| < 1e-8",
    "Variance Gamma max |phi(-i)-1| < 1e-8"
  ),
  value = c(
    as.numeric(isTRUE(black_sanity$summary$pass[1])),
    as.numeric(isTRUE(rough_cal$best$priced_eval$ok)),
    as.numeric(isTRUE(heston_cal$best$priced_eval$ok)),
    as.numeric(isTRUE(vg_cal$best$priced_eval$ok)),
    sum(rough_panel$model_price_clipped, na.rm = TRUE),
    sum(heston_panel$model_price_clipped, na.rm = TRUE),
    sum(vg_panel$model_price_clipped, na.rm = TRUE),
    safe_max(rheston_martingale_checks$abs_phi0_minus_1),
    safe_max(rheston_martingale_checks$abs_phim1_minus_1),
    safe_max(vg_martingale_checks$abs_phi0_minus_1),
    safe_max(vg_martingale_checks$abs_phim1_minus_1)
  ),
  pass = c(
    isTRUE(black_sanity$summary$pass[1]),
    isTRUE(rough_cal$best$priced_eval$ok),
    isTRUE(heston_cal$best$priced_eval$ok),
    isTRUE(vg_cal$best$priced_eval$ok),
    sum(rough_panel$model_price_clipped, na.rm = TRUE) == 0,
    sum(heston_panel$model_price_clipped, na.rm = TRUE) == 0,
    sum(vg_panel$model_price_clipped, na.rm = TRUE) == 0,
    safe_max(rheston_martingale_checks$abs_phi0_minus_1) < 1e-6,
    safe_max(rheston_martingale_checks$abs_phim1_minus_1) < 1e-6,
    safe_max(vg_martingale_checks$abs_phi0_minus_1) < 1e-8,
    safe_max(vg_martingale_checks$abs_phim1_minus_1) < 1e-8
  )
)

write_csv(acceptance, file.path(out_dir, "phaseM2010_empirical_calibration_acceptance_checks.csv"))
print(acceptance)
print(ranking_summary)

# ============================================================
# 9. Plots
# ============================================================

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 8.5),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

plot_days <- selected_expiries

fit_plot_long <- fit_panel %>%
  filter(tau_days %in% plot_days, sample %in% c("calibration", "holdout")) %>%
  select(tau_days, tau, k, iv_mkt, model, model_role, model_iv, sample) %>%
  mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(plot_days, "d")))

p_fit <- fit_plot_long %>%
  ggplot(aes(x = k)) +
  geom_point(aes(y = iv_mkt, shape = sample), size = 1.15, alpha = 0.75) +
  geom_line(aes(y = model_iv, colour = model, linetype = model_role, group = model), linewidth = 0.55) +
  facet_wrap(~ tau_label, scales = "free_y") +
  labs(
    x = "Forward log-moneyness, k = log(K/F)",
    y = "B&S implied volatility",
    colour = "Model",
    linetype = "Role",
    shape = "Sample"
  ) +
  plot_theme

plot_save(p_fit, "plot_phaseM2010_market_vs_calibrated_models", width = 10.5, height = 6.2)

p_resid <- fit_panel %>%
  filter(tau_days %in% plot_days, sample %in% c("calibration", "holdout")) %>%
  mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(plot_days, "d"))) %>%
  ggplot(aes(x = k, y = iv_error, colour = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_line(linewidth = 0.45) +
  geom_point(aes(shape = sample), size = 1.1) +
  facet_wrap(~ tau_label, scales = "free_y") +
  labs(
    x = "Forward log-moneyness, k = log(K/F)",
    y = "Model B&S-IV - market B&S-IV",
    colour = "Model",
    shape = "Sample"
  ) +
  plot_theme

plot_save(p_resid, "plot_phaseM2010_calibration_iv_residuals", width = 10.5, height = 6.2)

skew_plot <- bind_rows(
  market_skew_eval %>%
    select(tau_days, tau, atm_skew) %>%
    mutate(series = "market", model_role = "market"),
  atm_skew_comparison %>%
    select(tau_days, tau, atm_skew, model) %>%
    mutate(
      series = model,
      model_role = ifelse(model == "rough Heston", "main", "benchmark")
    ) %>%
    select(-model)
) %>%
  filter(is.finite(atm_skew), tau_days %in% selected_expiries)

p_skew <- skew_plot %>%
  ggplot(aes(x = tau_days, y = atm_skew, colour = series, linetype = model_role, shape = series)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  labs(
    x = "Maturity in calendar days",
    y = expression(paste("Signed ATM skew, ", Psi(tau))),
    colour = NULL,
    linetype = "Role",
    shape = NULL
  ) +
  plot_theme

plot_save(p_skew, "plot_phaseM2010_atm_skew_market_vs_models", width = 8.0, height = 5.0)

p_region <- error_by_region %>%
  filter(sample == "holdout") %>%
  ggplot(aes(x = moneyness_region, y = iv_rmse, fill = model)) +
  geom_col(position = "dodge", width = 0.72) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Holdout B&S-IV RMSE",
    fill = "Model"
  ) +
  plot_theme

plot_save(p_region, "plot_phaseM2010_holdout_iv_rmse_by_region", width = 8.2, height = 4.8)

figure_manifest <- tibble(
  file = paste0(c(
    "plot_phaseM2010_market_vs_calibrated_models",
    "plot_phaseM2010_calibration_iv_residuals",
    "plot_phaseM2010_atm_skew_market_vs_models",
    "plot_phaseM2010_holdout_iv_rmse_by_region"
  ), ".", PLOT_FORMAT),
  include_in_report = c(TRUE, FALSE, TRUE, TRUE),
  description = c(
    "Market B&S-IV versus calibrated rough Heston, classical Heston, and Variance Gamma models.",
    "Calibration B&S-IV residuals by maturity and model.",
    "Market signed ATM-skew term structure versus the main rough Heston model and benchmark models.",
    "Holdout B&S implied-volatility RMSE by moneyness region, showing model limitations away from ATM."
  )
)
write_csv(figure_manifest, file.path(out_dir, "phaseM2010_figure_manifest.csv"))

cat("\n=== Phase M2 2010 empirical model comparison completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("FAST_TEST mode:", FAST_TEST, "\n")
cat("EXTENDED_ROBUSTNESS mode:", EXTENDED_ROBUSTNESS, "\n")
cat("Lambda skew in objective:", LAMBDA_SKEW, "\n")
cat("Selected expiries:", paste(selected_expiries, collapse = ", "), "days\n")
cat("Calibration quotes:", nrow(calib_panel), "\n")
cat("Holdout quotes:", nrow(holdout_panel), "\n")
cat("\nBest parameters:\n")
print(best_params)
cat("\nError summary:\n")
print(error_summary)
cat("\nATM-skew summary:\n")
print(atm_skew_summary)
cat("\nModel ranking summary:\n")
print(ranking_summary)
cat("\nAcceptance checks:\n")
print(acceptance)
cat("\nMain files to inspect:\n")
cat("  phaseM2010_start_design_summary.csv\n")
cat("  phaseM2010_optimization_robustness_summary.csv\n")
cat("  phaseM2010_calibrated_parameters_best.csv\n")
cat("  phaseM2010_model_error_summary.csv\n")
cat("  phaseM2010_model_error_by_region.csv\n")
cat("  phaseM2010_atm_skew_market_vs_models.csv\n")
cat("  phaseM2010_atm_skew_error_summary.csv\n")
cat("  phaseM2010_atm_skew_error_summary_headline_H_maturities.csv\n")
cat("  phaseM2010_atm_skew_extraction_support_audit.csv\n")
cat("  phaseM2010_model_ranking_summary.csv\n")
cat("  phaseM2010_empirical_calibration_acceptance_checks.csv\n")
