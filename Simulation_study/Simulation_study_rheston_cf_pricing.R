#!/usr/bin/env Rscript

# ============================================================
# Phase S pricing: rough Heston characteristic-function validation
# ------------------------------------------------------------
# File:
#   Phase_s_rheston_cf_pricing.R
#
# Purpose:
#   Price the same Volterra Heston grid used in the iVi Monte Carlo
#   reference and compare the Fourier/Lewis output against that
#   independent simulation benchmark.
#
# Main diagnostics:
#   1. Black--Scholes/Lewis control check.
#   2. Characteristic-function normalization and martingale checks.
#   3. No-clipping audit in the near-ATM reporting region.
#   4. ATM-skew comparison between CF/Lewis and iVi.
#   5. alpha = 1 classical Heston boundary check.
#   6. Fractional-Riccati time-grid convergence.
#   7. Lewis/Fourier truncation convergence.
#   8. H-sensitivity of the ATM-skew term structure.
#
# Reduced run mode:
#   Set PHASES_FAST_TEST=1 after running the simulation script in
#   the same mode.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(161803)

# ============================================================
# 0. Run configuration and paths
# ============================================================

FAST_TEST <- identical(Sys.getenv("PHASES_FAST_TEST"), "1")

CLEAN_OUTPUT_DIR <- TRUE
PLOT_FORMAT <- "pdf"
SHOW_PLOTS <- !FAST_TEST
SAVE_PLOTS <- TRUE

CORE_K_MAX <- 0.05
STRESS_TAU_DAYS <- 7
STRESS_K <- 0.10

DU_TARGET <- 0.125

u_max_for_tau_days <- function(td) {
  # Fourier-grid guardrail: reduced mode may use fewer maturities and
  # Riccati steps, but the Lewis control check still needs adequate
  # resolution at short maturities.
  if (FAST_TEST) {
    if (td <= 7) return(180)
    if (td <= 30) return(140)
    return(120)
  }
  if (td <= 7) return(280)
  if (td <= 14) return(220)
  if (td <= 30) return(180)
  if (td <= 90) return(150)
  return(140)
}

n_time_for_tau_days <- function(td) {
  if (FAST_TEST) {
    if (td <= 7) return(80)
    if (td <= 30) return(70)
    return(60)
  }
  if (td <= 7) return(1100)
  if (td <= 14) return(800)
  if (td <= 30) return(620)
  return(520)
}

RICCATI_CORRECTOR_ITERS <- 2

RUN_ALPHA_ONE_LIMIT <- TRUE
RUN_RICCATI_TIME_CONVERGENCE <- !FAST_TEST
RUN_FOURIER_CONVERGENCE <- !FAST_TEST
RUN_H_SENSITIVITY <- !FAST_TEST

if (FAST_TEST) {
  RICCATI_CONV_N_TIME_GRID <- c(40, 60, 80)
  RICCATI_CONV_TAU_DAYS <- c(7, 30)
  FOURIER_CONV_U_MAX_GRID <- c(40, 60, 80)
  SENSITIVITY_H_GRID <- c(0.10, 0.50)
} else {
  RICCATI_CONV_N_TIME_GRID <- c(400, 620, 800, 1100, 1400)
  RICCATI_CONV_TAU_DAYS <- c(7, 30, 180)
  FOURIER_CONV_U_MAX_GRID <- c(100, 140, 180, 220, 280, 340)
  SENSITIVITY_H_GRID <- c(0.05, 0.10, 0.20, 0.50)
}

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

sim_dir <- file.path(base_dir, "data", "SimulationStudy", "phaseS_v2_output")
out_dir <- file.path(base_dir, "data", "SimulationStudy", "phaseS_rheston_cf_pricing_output")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_OUTPUT_DIR) {
  old_files <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)
}

cat("Input simulation directory:\n", sim_dir, "\n\n", sep = "")
cat("Output pricing directory:\n", out_dir, "\n\n", sep = "")
cat("FAST_TEST mode:", FAST_TEST, "\n\n")

required_files <- c(
  "simulation_v2_ivi_parameters.csv",
  "simulation_v2_ivi_mc_reference_surface.csv",
  "simulation_v2_ivi_mc_reference_diagnostics.csv"
)

missing_files <- required_files[!file.exists(file.path(sim_dir, required_files))]
if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing required Phase S output files:\n",
      paste(missing_files, collapse = "\n"),
      "\nRun Phase_s_simulation_study.R first."
    )
  )
}

plot_save <- function(plot_obj, stem, width, height, dpi = 180) {
  if (isTRUE(SHOW_PLOTS)) print(plot_obj)
  if (isTRUE(SAVE_PLOTS)) {
    path <- file.path(out_dir, paste0(stem, ".", PLOT_FORMAT))
    if (PLOT_FORMAT == "png") {
      ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi)
    } else if (PLOT_FORMAT == "pdf") {
      ggsave(filename = path, plot = plot_obj, width = width, height = height)
    } else {
      stop("PLOT_FORMAT must be either 'pdf' or 'png'.")
    }
  }
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
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
# 1. Read iVi reference outputs
# ============================================================

ivi_parameters <- read_csv(file.path(sim_dir, "simulation_v2_ivi_parameters.csv"), show_col_types = FALSE)
ivi_reference_surface <- read_csv(file.path(sim_dir, "simulation_v2_ivi_mc_reference_surface.csv"), show_col_types = FALSE)
ivi_reference_diagnostics <- read_csv(file.path(sim_dir, "simulation_v2_ivi_mc_reference_diagnostics.csv"), show_col_types = FALSE)

parameter_scenarios_path <- file.path(sim_dir, "simulation_v2_parameter_scenarios.csv")
parameter_scenarios <- if (file.exists(parameter_scenarios_path)) {
  read_csv(parameter_scenarios_path, show_col_types = FALSE)
} else {
  tibble()
}

F0 <- median(ivi_reference_surface$K / exp(ivi_reference_surface$k), na.rm = TRUE)
D0 <- 1

# Model parameters inherited from the iVi reference run.
rh_par <- list(
  gamma = ivi_parameters$gamma[1],
  theta = ivi_parameters$theta[1],
  nu = ivi_parameters$nu[1],
  rho = ivi_parameters$rho[1],
  V0 = ivi_parameters$V0[1],
  alpha = ivi_parameters$alpha[1]
)
rh_par$H <- rh_par$alpha - 0.5

pricing_grid <- ivi_reference_surface %>%
  distinct(tau_days, tau, k, K) %>%
  arrange(tau_days, k)

maturity_list <- sort(unique(pricing_grid$tau_days))

u_rule_tbl <- tibble(
  tau_days = maturity_list,
  u_max = sapply(maturity_list, u_max_for_tau_days),
  n_u = ceiling(u_max / DU_TARGET),
  du_target = DU_TARGET,
  riccati_n_time = sapply(maturity_list, n_time_for_tau_days)
)

parameter_tbl <- tibble(
  scenario_id = if ("scenario_id" %in% names(ivi_parameters)) ivi_parameters$scenario_id[1] else "iVi_stable_base",
  parameter_status = if ("parameter_status" %in% names(ivi_parameters)) ivi_parameters$parameter_status[1] else "model_based_simulation_benchmark_not_market_calibrated",
  gamma = rh_par$gamma,
  theta = rh_par$theta,
  nu = rh_par$nu,
  rho = rh_par$rho,
  V0 = rh_par$V0,
  alpha = rh_par$alpha,
  H = rh_par$H,
  F0 = F0,
  D0 = D0,
  du_target = DU_TARGET,
  riccati_corrector_iters = RICCATI_CORRECTOR_ITERS,
  core_k_max = CORE_K_MAX,
  fast_test = FAST_TEST
)

write_csv(parameter_tbl, file.path(out_dir, "phaseS_rheston_cf_pricing_parameters.csv"))
write_csv(u_rule_tbl, file.path(out_dir, "phaseS_rheston_cf_u_rule_by_maturity.csv"))
write_csv(pricing_grid, file.path(out_dir, "phaseS_rheston_cf_pricing_grid.csv"))
if (nrow(parameter_scenarios) > 0) write_csv(parameter_scenarios, file.path(out_dir, "phaseS_parameter_scenarios_inherited.csv"))

print(parameter_tbl)
print(u_rule_tbl)

# ============================================================
# 2. Black--Scholes functions in forward form
# ============================================================

# Prices are expressed as D_tau times a forward-price expectation.
black_call_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) return(D * max(F - K, 0))
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

implied_vol_black_call <- function(price, F, K, tau, D = 1, tol = 1e-9) {
  intrinsic <- D * max(F - K, 0)
  upper <- D * F
  if (!is.finite(price)) return(NA_real_)
  if (price < intrinsic - 1e-8 || price > upper + 1e-8) return(NA_real_)
  price <- min(max(price, intrinsic + 1e-12), upper - 1e-12)
  f <- function(sig) black_call_forward(F, K, tau, sig, D) - price
  flo <- f(1e-5)
  fhi <- f(5.0)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  uniroot(f, lower = 1e-5, upper = 5.0, tol = tol)$root
}

# ============================================================
# 3. Lewis inversion and Black--Scholes control check
# ============================================================

# Characteristic function of the log-forward return under Black--Scholes.
black_cf_log_forward <- function(z, tau, sigma) {
  exp(-0.5 * sigma^2 * tau * (z^2 + 1i * z))
}

# Lewis formula on the Im(z) = -1/2 shifted contour.
lewis_call_from_cf_generic <- function(k, cf_vec, u_grid, F = F0, D = D0) {
  K <- F * exp(k)
  du <- u_grid[2] - u_grid[1]
  w <- rep(du, length(u_grid))
  w[1L] <- 0.5 * du
  w[length(w)] <- 0.5 * du
  vals <- exp(-1i * u_grid * k) * cf_vec / (u_grid^2 + 0.25)
  integral <- sum(w * Re(vals))
  Re(D * (F - sqrt(F * K) / pi * integral))
}

run_black_lewis_sanity <- function() {
  test_grid <- expand.grid(
    k = c(-0.10, -0.05, 0, 0.05, 0.10),
    tau_days = maturity_list
  ) %>%
    as_tibble() %>%
    mutate(tau = tau_days / 365)
  
  sigma <- 0.22
  out <- vector("list", nrow(test_grid))
  for (i in seq_len(nrow(test_grid))) {
    k_i <- test_grid$k[i]
    tau_i <- test_grid$tau[i]
    td_i <- test_grid$tau_days[i]
    K_i <- F0 * exp(k_i)
    u_max_i <- u_max_for_tau_days(td_i)
    n_u_i <- ceiling(u_max_i / DU_TARGET)
    u_grid_i <- seq(0, u_max_i, length.out = n_u_i + 1L)
    cf_vec <- sapply(u_grid_i, function(u) black_cf_log_forward(u - 0.5i, tau_i, sigma))
    p_black <- black_call_forward(F0, K_i, tau_i, sigma, D0)
    p_lewis <- lewis_call_from_cf_generic(k_i, cf_vec, u_grid_i, F0, D0)
    out[[i]] <- tibble(k = k_i, tau = tau_i, tau_days = td_i, u_max = u_max_i, n_u = n_u_i, price_black = p_black, price_lewis = p_lewis, abs_error = abs(p_black - p_lewis))
  }
  details <- bind_rows(out)
  summary <- details %>%
    summarise(max_abs_error = max(abs_error, na.rm = TRUE), mean_abs_error = mean(abs_error, na.rm = TRUE), pass = max_abs_error < 1e-4)
  list(details = details, summary = summary)
}

black_sanity <- run_black_lewis_sanity()
write_csv(black_sanity$details, file.path(out_dir, "phaseS_black_lewis_sanity_details.csv"))
write_csv(black_sanity$summary, file.path(out_dir, "phaseS_black_lewis_sanity.csv"))
print(black_sanity$summary)
if (!isTRUE(black_sanity$summary$pass[1])) {
  print(black_sanity$details %>% arrange(desc(abs_error)) %>% head(10))
  stop("Black/Lewis sanity check failed. Fix Lewis grid before rough-Heston pricing.")
}

# ============================================================
# 4. Fractional Riccati solver and alpha = 1 boundary case
# ============================================================

is_finite_complex <- function(z) is.finite(Re(z)) && is.finite(Im(z))

# Adams--Bashforth--Moulton weights for the fractional Riccati equation.
make_abm_precomp <- function(T, n_time, alpha) {
  if (alpha >= 0.999999) stop("ABM precomputation is for alpha < 1. Use classical Heston ODE for alpha=1.")
  dt <- T / n_time
  n <- n_time
  lag <- 1:n
  b_lag <- (dt^alpha / base::gamma(alpha + 1)) * (lag^alpha - (lag - 1)^alpha)
  a_list <- vector("list", n)
  for (k in 1:n) {
    a <- numeric(k)
    a[1] <- (k - 1)^(alpha + 1) - (k - alpha - 1) * k^alpha
    if (k >= 2) {
      for (j in 1:(k - 1)) {
        lag_kj <- k - j
        a[j + 1] <- (lag_kj + 1)^(alpha + 1) - 2 * lag_kj^(alpha + 1) + (lag_kj - 1)^(alpha + 1)
      }
    }
    a_list[[k]] <- (dt^alpha / base::gamma(alpha + 2)) * a
  }
  beta <- 1 - alpha
  lag_beta <- n:1
  weights_beta <- ((lag_beta^beta - (lag_beta - 1)^beta) * dt^beta) / base::gamma(beta + 1)
  list(T = T, n_time = n_time, alpha = alpha, dt = dt, b_lag = b_lag, a_list = a_list, weights_beta = weights_beta)
}

riccati_rhs <- function(a, y, par) {
  gamma_p <- par$gamma
  nu_p <- par$nu
  rho_p <- par$rho
  ccoef <- gamma_p * nu_p
  0.5 * (-a^2 - 1i * a) +
    gamma_p * (1i * a * rho_p * nu_p - 1) * y +
    0.5 * ccoef^2 * y^2
}

# Rough Heston transform for one complex argument.
rough_heston_cf_one_abm <- function(a, pre, par, corrector_iters = 2) {
  n <- pre$n_time
  dt <- pre$dt
  gamma_p <- par$gamma
  theta_p <- par$theta
  V0_p <- par$V0
  
  h <- complex(length = n + 1L)
  fvals <- complex(length = n + 1L)
  h[1L] <- 0 + 0i
  fvals[1L] <- riccati_rhs(a, h[1L], par)
  
  ok <- TRUE
  fail_step <- NA_integer_
  for (k in 1:n) {
    h_pred <- sum(pre$b_lag[k:1] * fvals[seq_len(k)])
    if (!is_finite_complex(h_pred) || Mod(h_pred) > 1e10) {
      ok <- FALSE; fail_step <- k; break
    }
    a_weights <- pre$a_list[[k]]
    base_term <- sum(a_weights * fvals[seq_len(k)])
    h_new <- h_pred
    for (it in seq_len(corrector_iters)) {
      f_new <- riccati_rhs(a, h_new, par)
      h_new <- base_term + (dt^pre$alpha / base::gamma(pre$alpha + 2)) * f_new
      if (!is_finite_complex(h_new) || Mod(h_new) > 1e10) {
        ok <- FALSE; fail_step <- k; break
      }
    }
    if (!ok) break
    h[k + 1L] <- h_new
    fvals[k + 1L] <- riccati_rhs(a, h_new, par)
    if (!is_finite_complex(fvals[k + 1L]) || Mod(fvals[k + 1L]) > 1e12) {
      ok <- FALSE; fail_step <- k; break
    }
  }
  
  if (!ok) {
    return(list(cf = NA_complex_, exponent = NA_complex_, ok = FALSE, fail_step = fail_step, max_mod_h = safe_max(Mod(h))))
  }
  
  int_h <- dt * (sum(h) - 0.5 * h[1L] - 0.5 * h[length(h)])
  I_beta_h_T <- sum(pre$weights_beta * h[seq_len(n)])
  exponent <- gamma_p * theta_p * int_h + V0_p * I_beta_h_T
  if (!is_finite_complex(exponent) || Re(exponent) > 700) {
    return(list(cf = NA_complex_, exponent = exponent, ok = FALSE, fail_step = NA_integer_, max_mod_h = safe_max(Mod(h))))
  }
  cf <- exp(exponent)
  list(cf = cf, exponent = exponent, ok = is_finite_complex(cf), fail_step = NA_integer_, max_mod_h = safe_max(Mod(h)))
}

# Classical Heston boundary case obtained when alpha = 1.
classical_heston_cf_one_ode <- function(a, T, par, n_time = 2000) {
  n <- max(10L, as.integer(n_time))
  dt <- T / n
  h <- 0 + 0i
  int_h <- 0 + 0i
  max_mod_h <- 0
  
  f <- function(y) riccati_rhs(a, y, par)
  
  for (i in seq_len(n)) {
    h_old <- h
    k1 <- f(h)
    k2 <- f(h + 0.5 * dt * k1)
    k3 <- f(h + 0.5 * dt * k2)
    k4 <- f(h + dt * k3)
    h <- h + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    if (!is_finite_complex(h) || Mod(h) > 1e10) {
      return(list(cf = NA_complex_, exponent = NA_complex_, ok = FALSE, fail_step = i, max_mod_h = max_mod_h))
    }
    int_h <- int_h + 0.5 * dt * (h_old + h)
    max_mod_h <- max(max_mod_h, Mod(h))
  }
  
  exponent <- par$gamma * par$theta * int_h + par$V0 * h
  cf <- exp(exponent)
  list(cf = cf, exponent = exponent, ok = is_finite_complex(cf), fail_step = NA_integer_, max_mod_h = max_mod_h)
}

rheston_cf_one <- function(a, T, par, n_time, corrector_iters = 2) {
  if (par$alpha >= 0.999999) {
    classical_heston_cf_one_ode(a, T, par, n_time = n_time)
  } else {
    pre <- make_abm_precomp(T = T, n_time = n_time, alpha = par$alpha)
    rough_heston_cf_one_abm(a = a, pre = pre, par = par, corrector_iters = corrector_iters)
  }
}

compute_cf_grid_for_maturity <- function(T, u_grid, par, n_time, corrector_iters = 2, progress = TRUE) {
  pre <- NULL
  if (par$alpha < 0.999999) {
    pre <- make_abm_precomp(T = T, n_time = n_time, alpha = par$alpha)
  }
  cf_out <- vector("list", length(u_grid))
  for (idx in seq_along(u_grid)) {
    if (progress && idx %% 300 == 0) cat("      CF u-index", idx, "of", length(u_grid), "\n")
    a <- u_grid[idx] - 0.5i
    res <- if (par$alpha >= 0.999999) {
      classical_heston_cf_one_ode(a = a, T = T, par = par, n_time = n_time)
    } else {
      rough_heston_cf_one_abm(a = a, pre = pre, par = par, corrector_iters = corrector_iters)
    }
    cf_out[[idx]] <- tibble(u = u_grid[idx], cf_real = Re(res$cf), cf_imag = Im(res$cf), cf_ok = isTRUE(res$ok), max_mod_h = res$max_mod_h, fail_step = res$fail_step)
  }
  bind_rows(cf_out)
}

# ============================================================
# 5. CF normalization / martingale checks
# ============================================================

martingale_checks <- lapply(maturity_list, function(td) {
  Tval <- td / 365
  nt <- n_time_for_tau_days(td)
  cf0 <- rheston_cf_one(a = 0 + 0i, T = Tval, par = rh_par, n_time = nt, corrector_iters = RICCATI_CORRECTOR_ITERS)
  cf_minus_i <- rheston_cf_one(a = -1i, T = Tval, par = rh_par, n_time = nt, corrector_iters = RICCATI_CORRECTOR_ITERS)
  tibble(
    tau_days = td,
    tau = Tval,
    n_time = nt,
    cf0_real = Re(cf0$cf),
    cf0_imag = Im(cf0$cf),
    martingale_cf_minus_i_real = Re(cf_minus_i$cf),
    martingale_cf_minus_i_imag = Im(cf_minus_i$cf),
    abs_cf0_minus_1 = Mod(cf0$cf - 1),
    abs_cf_minus_i_minus_1 = Mod(cf_minus_i$cf - 1),
    cf0_ok = cf0$ok,
    cf_minus_i_ok = cf_minus_i$ok
  )
}) %>% bind_rows()

write_csv(martingale_checks, file.path(out_dir, "phaseS_rheston_cf_martingale_checks.csv"))
print(martingale_checks)

if (max(martingale_checks$abs_cf0_minus_1, na.rm = TRUE) > 1e-8 ||
    max(martingale_checks$abs_cf_minus_i_minus_1, na.rm = TRUE) > 1e-8) {
  warning("CF normalization or martingale check is not tight.")
}

# ============================================================
# 6. Lewis pricing from the rough Heston characteristic function
# ============================================================

price_from_cf_for_k <- function(k, tau, cf_vec, u_grid, F = F0, D = D0) {
  K <- F * exp(k)
  price_raw <- lewis_call_from_cf_generic(k = k, cf_vec = cf_vec, u_grid = u_grid, F = F, D = D)
  intrinsic <- D * max(F - K, 0)
  upper <- D * F
  price_clipped <- min(max(price_raw, intrinsic + 1e-10), upper - 1e-10)
  clipped <- abs(price_raw - price_clipped) > 1e-7
  iv <- implied_vol_black_call(price_clipped, F, K, tau, D)
  vega <- ifelse(is.finite(iv), black_vega_forward(F, K, tau, iv, D), NA_real_)
  tibble(
    k = k,
    K = K,
    cf_lewis_price_raw = price_raw,
    cf_lewis_price = price_clipped,
    cf_price_clipped = clipped,
    cf_lewis_iv = iv,
    cf_black_vega = vega,
    intrinsic = intrinsic,
    upper = upper
  )
}

# Maturity-by-maturity pricing keeps the Fourier and Riccati grids maturity dependent.
price_surface_from_par <- function(par, grid, du_target = DU_TARGET, u_max_fun = u_max_for_tau_days, n_time_fun = n_time_for_tau_days, k_filter = NULL, progress = TRUE) {
  all_cf_tables <- list()
  all_price_tables <- list()
  maturity_vec <- sort(unique(grid$tau_days))
  for (td in maturity_vec) {
    Tval <- unique(grid$tau[grid$tau_days == td])[1]
    u_max_td <- u_max_fun(td)
    n_u_td <- ceiling(u_max_td / du_target)
    n_time_td <- n_time_fun(td)
    u_grid_td <- seq(0, u_max_td, length.out = n_u_td + 1L)
    
    if (progress) {
      cat("Pricing maturity tau_days =", td, "| Umax =", u_max_td, "| N_u =", n_u_td, "| Riccati N_time =", n_time_td, "\n")
    }
    
    cf_tbl <- compute_cf_grid_for_maturity(
      T = Tval,
      u_grid = u_grid_td,
      par = par,
      n_time = n_time_td,
      corrector_iters = RICCATI_CORRECTOR_ITERS,
      progress = progress
    ) %>%
      mutate(tau_days = td, tau = Tval, n_time = n_time_td, corrector_iters = RICCATI_CORRECTOR_ITERS, u_max = u_max_td, n_u = n_u_td)
    
    all_cf_tables[[as.character(td)]] <- cf_tbl
    cf_vec <- cf_tbl$cf_real + 1i * cf_tbl$cf_imag
    
    k_vec <- grid %>% filter(tau_days == td) %>% arrange(k) %>% pull(k)
    if (!is.null(k_filter)) k_vec <- k_vec[k_vec %in% k_filter]
    
    price_tbl <- bind_rows(lapply(k_vec, function(kval) {
      price_from_cf_for_k(k = kval, tau = Tval, cf_vec = cf_vec, u_grid = u_grid_td, F = F0, D = D0)
    })) %>%
      mutate(tau_days = td, tau = Tval, n_time = n_time_td, corrector_iters = RICCATI_CORRECTOR_ITERS, u_max = u_max_td, n_u = n_u_td) %>%
      select(tau_days, tau, k, K, everything())
    
    all_price_tables[[as.character(td)]] <- price_tbl
  }
  list(cf = bind_rows(all_cf_tables), prices = bind_rows(all_price_tables))
}

main_pricing <- price_surface_from_par(rh_par, pricing_grid, progress = TRUE)
cf_eval_grid <- main_pricing$cf
rheston_cf_prices <- main_pricing$prices

write_csv(cf_eval_grid, file.path(out_dir, "phaseS_rheston_cf_evaluation_grid.csv"))
write_csv(rheston_cf_prices, file.path(out_dir, "phaseS_rheston_cf_lewis_prices.csv"))

cf_grid_summary <- cf_eval_grid %>%
  group_by(tau_days) %>%
  summarise(
    n_u = n(),
    u_max = max(u_max),
    n_time = max(n_time),
    n_bad_cf = sum(!cf_ok | !is.finite(cf_real) | !is.finite(cf_imag)),
    first_bad_u = {
      bad_u <- u[!cf_ok | !is.finite(cf_real) | !is.finite(cf_imag)]
      if (length(bad_u) == 0L) NA_real_ else min(bad_u)
    },
    max_good_u = {
      good_u <- u[cf_ok & is.finite(cf_real) & is.finite(cf_imag)]
      if (length(good_u) == 0L) NA_real_ else max(good_u)
    },
    max_mod_h = safe_max(max_mod_h),
    .groups = "drop"
  )

write_csv(cf_grid_summary, file.path(out_dir, "phaseS_rheston_cf_grid_summary.csv"))
print(cf_grid_summary)

# ============================================================
# 7. CF/Lewis comparison against the iVi reference
# ============================================================

cf_vs_ivi <- rheston_cf_prices %>%
  left_join(
    ivi_reference_surface %>%
      select(
        tau_days, k,
        ivi_ref_tau = tau,
        ivi_ref_K = K,
        ivi_n_steps = n_steps,
        ivi_n_paths = n_paths,
        ivi_ref_price = mc_price,
        ivi_ref_price_se = mc_price_se,
        ivi_ref_iv = mc_iv,
        ivi_ref_iv_se = mc_iv_se
      ),
    by = c("tau_days", "k")
  ) %>%
  mutate(
    K_difference_cf_minus_ivi = K - ivi_ref_K,
    price_error_cf_minus_ivi = cf_lewis_price - ivi_ref_price,
    abs_price_error = abs(price_error_cf_minus_ivi),
    price_error_over_mcse = price_error_cf_minus_ivi / ivi_ref_price_se,
    iv_error_cf_minus_ivi = cf_lewis_iv - ivi_ref_iv,
    abs_iv_error = abs(iv_error_cf_minus_ivi),
    iv_error_over_mcse = iv_error_cf_minus_ivi / ivi_ref_iv_se,
    abs_iv_error_over_mcse = abs(iv_error_over_mcse),
    is_core_k = abs(k) <= CORE_K_MAX + 1e-12,
    is_stress_7d_right_wing = tau_days == STRESS_TAU_DAYS & abs(k - STRESS_K) <= 1e-12,
    validation_region = case_when(
      is_core_k ~ "headline_core_abs_k_le_0p05",
      is_stress_7d_right_wing ~ "stress_7d_right_wing_k_0p10",
      TRUE ~ "outer_wing_nonheadline"
    )
  )

write_csv(cf_vs_ivi, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_reference.csv"))

core_clip_count <- cf_vs_ivi %>%
  filter(is_core_k) %>%
  summarise(n_core_clipped = sum(cf_price_clipped, na.rm = TRUE)) %>%
  pull(n_core_clipped)
if (core_clip_count > 0L) warning("CF price clipping occurred inside the headline ATM/near-ATM region.")

# ============================================================
# 8. Summary tables
# ============================================================

summarise_errors <- function(df) {
  df %>%
    summarise(
      n_quotes = n(),
      n_cf_price_clipped = sum(cf_price_clipped, na.rm = TRUE),
      n_na_iv = sum(!is.finite(cf_lewis_iv)),
      price_rmse = safe_rmse(price_error_cf_minus_ivi),
      price_mae = safe_mae(price_error_cf_minus_ivi),
      mean_ivi_price_se = safe_mean(ivi_ref_price_se),
      price_rmse_over_mean_mcse = price_rmse / mean_ivi_price_se,
      iv_rmse = safe_rmse(iv_error_cf_minus_ivi),
      iv_mae = safe_mae(iv_error_cf_minus_ivi),
      mean_ivi_iv_se = safe_mean(ivi_ref_iv_se),
      iv_rmse_over_mean_mcse = iv_rmse / mean_ivi_iv_se,
      max_abs_price_error = safe_max(abs_price_error),
      max_abs_price_error_over_mcse = safe_max(abs(price_error_over_mcse)),
      max_abs_iv_error = safe_max(abs_iv_error),
      max_abs_iv_error_over_mcse = safe_max(abs_iv_error_over_mcse),
      .groups = "drop"
    )
}

summary_by_maturity <- cf_vs_ivi %>% group_by(tau_days, tau) %>% summarise_errors() %>% ungroup()
summary_by_region <- cf_vs_ivi %>% group_by(validation_region) %>% summarise_errors() %>% ungroup()
summary_by_maturity_region <- cf_vs_ivi %>% group_by(tau_days, tau, validation_region) %>% summarise_errors() %>% ungroup()
summary_global_full <- cf_vs_ivi %>% summarise_errors() %>% mutate(validation_region = "full_grid_all_points") %>% select(validation_region, everything())
summary_global_core <- cf_vs_ivi %>% filter(is_core_k) %>% summarise_errors() %>% mutate(validation_region = "headline_core_abs_k_le_0p05") %>% select(validation_region, everything())
summary_global_stress <- cf_vs_ivi %>% filter(is_stress_7d_right_wing) %>% summarise_errors() %>% mutate(validation_region = "stress_7d_right_wing_k_0p10") %>% select(validation_region, everything())
summary_global_outer <- cf_vs_ivi %>% filter(validation_region == "outer_wing_nonheadline") %>% summarise_errors() %>% mutate(validation_region = "outer_wing_nonheadline") %>% select(validation_region, everything())

summary_global_report <- bind_rows(summary_global_core, summary_global_full, summary_global_stress, summary_global_outer)

write_csv(summary_by_maturity, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_summary_by_maturity.csv"))
write_csv(summary_by_region, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_summary_by_region.csv"))
write_csv(summary_by_maturity_region, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_summary_by_maturity_region.csv"))
write_csv(summary_global_report, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_summary_global_report.csv"))
write_csv(summary_global_full, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_summary_global.csv"))

print(summary_by_maturity)
print(summary_by_region)
print(summary_global_report)

# ============================================================
# 9. ATM-skew comparison between CF/Lewis and iVi
# ============================================================

estimate_atm_skew <- function(df, iv_col, window = 0.055) {
  # ATM skew is estimated as the local derivative of the IV slice at k = 0.
  # The headline estimator is two-sided local quadratic; the local-linear
  # fallback is kept only for reduced runs and audit tables.
  d <- df %>% filter(abs(k) <= window + 1e-12, is.finite(.data[[iv_col]]))
  if (nrow(d) < 3L) {
    return(tibble(
      psi = NA_real_,
      sigma_atm = NA_real_,
      n_local = nrow(d),
      fit_type = "insufficient"
    ))
  }
  d <- d %>% mutate(local_weight = pmax(0, 1 - abs(k) / window))
  if (nrow(d) >= 5L) {
    fit <- lm(as.formula(paste0(iv_col, " ~ k + I(k^2)")), data = d, weights = local_weight)
    cf <- coef(fit)
    return(tibble(
      psi = unname(cf["k"]),
      sigma_atm = unname(cf["(Intercept)"]),
      n_local = nrow(d),
      fit_type = "local_quadratic"
    ))
  }
  fit <- lm(as.formula(paste0(iv_col, " ~ k")), data = d, weights = local_weight)
  cf <- coef(fit)
  tibble(
    psi = unname(cf["k"]),
    sigma_atm = unname(cf["(Intercept)"]),
    n_local = nrow(d),
    fit_type = "local_linear"
  )
}

# Records whether the local ATM fit has observations on both sides of k = 0.
near_atm_support_audit <- function(df, window = 0.055) {
  d <- df %>% filter(abs(k) <= window + 1e-12)
  if (nrow(d) == 0L) {
    return(tibble(
      n_near_atm_points = 0L,
      k_support = NA_character_,
      has_left = FALSE,
      has_atm = FALSE,
      has_right = FALSE,
      two_sided_support = FALSE
    ))
  }
  tibble(
    n_near_atm_points = dplyr::n_distinct(d$k),
    k_support = paste(sprintf("%.3f", sort(unique(d$k))), collapse = ", "),
    has_left = any(d$k < -1e-12),
    has_atm = any(abs(d$k) <= 1e-12),
    has_right = any(d$k > 1e-12),
    two_sided_support = has_left & has_atm & has_right
  )
}

central_difference_skew <- function(df, iv_col, h = 0.025, tol = 1e-10) {
  d <- df %>% filter(is.finite(k), is.finite(.data[[iv_col]]))
  idx_plus <- which(abs(d$k - h) <= tol)
  idx_minus <- which(abs(d$k + h) <= tol)
  if (length(idx_plus) < 1L || length(idx_minus) < 1L) {
    return(tibble(cd_skew = NA_real_, cd_h = h, cd_available = FALSE))
  }
  iv_plus <- d[[iv_col]][idx_plus[1L]]
  iv_minus <- d[[iv_col]][idx_minus[1L]]
  tibble(
    cd_skew = (iv_plus - iv_minus) / (2 * h),
    cd_h = h,
    cd_available = TRUE
  )
}

# Recovers H from the slope in log |psi(tau)| against log tau.
fit_loglog_skew <- function(df, skew_col, method_label, H_target) {
  d <- df %>%
    filter(
      is.finite(.data[[skew_col]]),
      abs(.data[[skew_col]]) > 0,
      is.finite(tau),
      tau > 0
    ) %>%
    mutate(abs_skew = abs(.data[[skew_col]]))
  if (nrow(d) < 3L) {
    return(list(
      summary = tibble(
        method = method_label,
        intercept = NA_real_,
        beta_hat = NA_real_,
        H_recovered = NA_real_,
        target_H = H_target,
        target_beta = H_target - 0.5,
        r_squared = NA_real_,
        n_maturities = nrow(d)
      ),
      fitted = tibble()
    ))
  }
  fit <- lm(log(abs_skew) ~ log(tau), data = d)
  sfit <- summary(fit)
  fitted_tbl <- d %>%
    mutate(
      method = method_label,
      log_abs_skew = log(abs_skew),
      fitted_log_abs_skew = fitted(fit),
      loglog_residual = residuals(fit)
    ) %>%
    select(method, tau_days, tau, abs_skew, log_abs_skew, fitted_log_abs_skew, loglog_residual)
  summary_tbl <- tibble(
    method = method_label,
    intercept = unname(coef(fit)[1]),
    beta_hat = unname(coef(fit)[2]),
    H_recovered = beta_hat + 0.5,
    target_H = H_target,
    target_beta = H_target - 0.5,
    r_squared = sfit$r.squared,
    n_maturities = nrow(d)
  )
  list(summary = summary_tbl, fitted = fitted_tbl)
}

atm_cf <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "cf_lewis_iv", window = 0.055)) %>%
  ungroup() %>%
  rename(cf_atm_skew = psi, cf_atm_iv = sigma_atm, cf_n_local = n_local, cf_fit_type = fit_type)

atm_ivi <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "ivi_ref_iv", window = 0.055)) %>%
  ungroup() %>%
  rename(ivi_atm_skew = psi, ivi_atm_iv = sigma_atm, ivi_n_local = n_local, ivi_fit_type = fit_type)

atm_skew_comparison <- atm_cf %>%
  left_join(atm_ivi, by = c("tau_days", "tau")) %>%
  mutate(
    atm_iv_error_cf_minus_ivi = cf_atm_iv - ivi_atm_iv,
    atm_skew_error_cf_minus_ivi = cf_atm_skew - ivi_atm_skew,
    abs_atm_skew_error = abs(atm_skew_error_cf_minus_ivi)
  )

near_atm_support <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ near_atm_support_audit(.x, window = 0.055)) %>%
  ungroup()

atm_skew_extraction_audit <- atm_skew_comparison %>%
  left_join(near_atm_support, by = c("tau_days", "tau")) %>%
  mutate(
    cf_uses_headline_estimator = cf_fit_type == "local_quadratic" & cf_n_local >= 5L,
    ivi_uses_headline_estimator = ivi_fit_type == "local_quadratic" & ivi_n_local >= 5L,
    headline_skew_estimator_ok = two_sided_support & cf_uses_headline_estimator & ivi_uses_headline_estimator
  )

cd_cf_h0025 <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "cf_lewis_iv", h = 0.025)) %>%
  ungroup() %>%
  rename(cf_cd_skew_h0025 = cd_skew, cf_cd_available_h0025 = cd_available)

cd_ivi_h0025 <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "ivi_ref_iv", h = 0.025)) %>%
  ungroup() %>%
  rename(ivi_cd_skew_h0025 = cd_skew, ivi_cd_available_h0025 = cd_available)

cd_cf_h005 <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "cf_lewis_iv", h = 0.05)) %>%
  ungroup() %>%
  rename(cf_cd_skew_h005 = cd_skew, cf_cd_available_h005 = cd_available)

cd_ivi_h005 <- cf_vs_ivi %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "ivi_ref_iv", h = 0.05)) %>%
  ungroup() %>%
  rename(ivi_cd_skew_h005 = cd_skew, ivi_cd_available_h005 = cd_available)

atm_skew_central_difference_check <- atm_skew_comparison %>%
  left_join(cd_cf_h0025, by = c("tau_days", "tau")) %>%
  left_join(cd_ivi_h0025, by = c("tau_days", "tau")) %>%
  left_join(cd_cf_h005, by = c("tau_days", "tau")) %>%
  left_join(cd_ivi_h005, by = c("tau_days", "tau")) %>%
  mutate(
    cf_lq_minus_cd_h0025 = cf_atm_skew - cf_cd_skew_h0025,
    ivi_lq_minus_cd_h0025 = ivi_atm_skew - ivi_cd_skew_h0025,
    cf_lq_minus_cd_h005 = cf_atm_skew - cf_cd_skew_h005,
    ivi_lq_minus_cd_h005 = ivi_atm_skew - ivi_cd_skew_h005
  )

cf_loglog <- fit_loglog_skew(atm_skew_comparison, "cf_atm_skew", "CF/Lewis", rh_par$H)
ivi_loglog <- fit_loglog_skew(atm_skew_comparison, "ivi_atm_skew", "iVi reference", rh_par$H)
atm_skew_loglog_summary <- bind_rows(cf_loglog$summary, ivi_loglog$summary)
atm_skew_loglog_points <- bind_rows(cf_loglog$fitted, ivi_loglog$fitted)

write_csv(atm_skew_comparison, file.path(out_dir, "phaseS_rheston_cf_vs_ivi_atm_skew_comparison.csv"))
write_csv(atm_skew_extraction_audit, file.path(out_dir, "phaseS_rheston_atm_skew_extraction_audit.csv"))
write_csv(atm_skew_central_difference_check, file.path(out_dir, "phaseS_rheston_atm_skew_central_difference_check.csv"))
write_csv(atm_skew_loglog_summary, file.path(out_dir, "phaseS_rheston_atm_skew_loglog_summary.csv"))
write_csv(atm_skew_loglog_points, file.path(out_dir, "phaseS_rheston_atm_skew_loglog_points.csv"))
print(atm_skew_comparison)
print(atm_skew_extraction_audit)
print(atm_skew_central_difference_check)
print(atm_skew_loglog_summary)

if (!FAST_TEST && !isTRUE(all(atm_skew_extraction_audit$headline_skew_estimator_ok))) {
  stop("ATM-skew extraction audit failed: headline maturities do not all use two-sided local-quadratic support.")
}

# ============================================================
# 10. alpha=1 classical-Heston limit check
# ============================================================

if (RUN_ALPHA_ONE_LIMIT) {
  cat("Running alpha=1 classical-Heston limit check...\n")
  alpha_grid <- if (FAST_TEST) c(0.90, 1.00) else c(0.90, 0.95, 0.98, 0.995, 1.00)
  alpha_limit_tau_days <- intersect(c(7, 30, 180), maturity_list)
  if (length(alpha_limit_tau_days) == 0L) alpha_limit_tau_days <- maturity_list[1]
  rows <- list()
  idx <- 1L
  for (td in alpha_limit_tau_days) {
    Tval <- td / 365
    u_max_td <- u_max_for_tau_days(td)
    n_u_td <- ceiling(u_max_td / DU_TARGET)
    u_grid_td <- seq(0, u_max_td, length.out = n_u_td + 1L)
    nt_td <- n_time_for_tau_days(td)
    for (alpha_val in alpha_grid) {
      par_a <- rh_par
      par_a$alpha <- alpha_val
      par_a$H <- alpha_val - 0.5
      cf_tbl <- compute_cf_grid_for_maturity(Tval, u_grid_td, par_a, nt_td, RICCATI_CORRECTOR_ITERS, progress = FALSE)
      cf_vec <- cf_tbl$cf_real + 1i * cf_tbl$cf_imag
      price_tbl <- price_from_cf_for_k(k = 0.0, tau = Tval, cf_vec = cf_vec, u_grid = u_grid_td, F = F0, D = D0)
      rows[[idx]] <- price_tbl %>% mutate(tau_days = td, tau = Tval, alpha = alpha_val, H = alpha_val - 0.5, n_time = nt_td, u_max = u_max_td, n_u = n_u_td)
      idx <- idx + 1L
    }
  }
  alpha_one_limit <- bind_rows(rows) %>%
    group_by(tau_days) %>%
    mutate(
      reference_alpha = 1.0,
      reference_iv_alpha_1 = cf_lewis_iv[which.min(abs(alpha - 1.0))][1],
      reference_price_alpha_1 = cf_lewis_price[which.min(abs(alpha - 1.0))][1],
      iv_diff_to_alpha_1 = cf_lewis_iv - reference_iv_alpha_1,
      price_diff_to_alpha_1 = cf_lewis_price - reference_price_alpha_1
    ) %>%
    ungroup()
  write_csv(alpha_one_limit, file.path(out_dir, "phaseS_alpha_one_heston_limit_check.csv"))
  print(alpha_one_limit)
}

# ============================================================
# 11. Riccati time-grid convergence
# ============================================================

if (RUN_RICCATI_TIME_CONVERGENCE) {
  cat("Running fractional-Riccati time-grid convergence check...\n")
  conv_rows <- list()
  idx <- 1L
  conv_tau <- intersect(RICCATI_CONV_TAU_DAYS, maturity_list)
  conv_k <- intersect(c(-0.05, 0, 0.05), unique(pricing_grid$k))
  if (length(conv_k) == 0L) conv_k <- 0
  for (td in conv_tau) {
    Tval <- td / 365
    u_max_td <- u_max_for_tau_days(td)
    n_u_td <- ceiling(u_max_td / DU_TARGET)
    u_grid_td <- seq(0, u_max_td, length.out = n_u_td + 1L)
    for (nt in RICCATI_CONV_N_TIME_GRID) {
      cat("  Riccati convergence: tau_days =", td, "| n_time =", nt, "\n")
      cf_tbl <- compute_cf_grid_for_maturity(Tval, u_grid_td, rh_par, nt, RICCATI_CORRECTOR_ITERS, progress = FALSE)
      cf_vec <- cf_tbl$cf_real + 1i * cf_tbl$cf_imag
      for (kval in conv_k) {
        conv_rows[[idx]] <- price_from_cf_for_k(kval, Tval, cf_vec, u_grid_td, F0, D0) %>%
          mutate(tau_days = td, tau = Tval, n_time = nt, u_max = u_max_td, n_u = n_u_td)
        idx <- idx + 1L
      }
    }
  }
  riccati_conv <- bind_rows(conv_rows) %>%
    group_by(tau_days, k) %>%
    mutate(
      reference_n_time = max(n_time),
      reference_iv = cf_lewis_iv[n_time == max(n_time)][1],
      reference_price = cf_lewis_price[n_time == max(n_time)][1],
      iv_error_to_fine = cf_lewis_iv - reference_iv,
      price_error_to_fine = cf_lewis_price - reference_price
    ) %>%
    ungroup()
  write_csv(riccati_conv, file.path(out_dir, "phaseS_rheston_riccati_time_convergence.csv"))
  print(riccati_conv)
}

# ============================================================
# 12. Fourier truncation convergence
# ============================================================

if (RUN_FOURIER_CONVERGENCE) {
  cat("Running Lewis/Fourier truncation convergence check...\n")
  fourier_rows <- list()
  idx <- 1L
  conv_tau <- intersect(c(7, 30, 180), maturity_list)
  conv_k <- intersect(c(-0.05, 0, 0.05), unique(pricing_grid$k))
  if (length(conv_k) == 0L) conv_k <- 0
  for (td in conv_tau) {
    Tval <- td / 365
    nt_td <- n_time_for_tau_days(td)
    for (u_max_val in FOURIER_CONV_U_MAX_GRID) {
      cat("  Fourier convergence: tau_days =", td, "| Umax =", u_max_val, "\n")
      n_u_val <- ceiling(u_max_val / DU_TARGET)
      u_grid_val <- seq(0, u_max_val, length.out = n_u_val + 1L)
      cf_tbl <- compute_cf_grid_for_maturity(Tval, u_grid_val, rh_par, nt_td, RICCATI_CORRECTOR_ITERS, progress = FALSE)
      cf_vec <- cf_tbl$cf_real + 1i * cf_tbl$cf_imag
      for (kval in conv_k) {
        fourier_rows[[idx]] <- price_from_cf_for_k(kval, Tval, cf_vec, u_grid_val, F0, D0) %>%
          mutate(tau_days = td, tau = Tval, u_max = u_max_val, n_u = n_u_val, n_time = nt_td)
        idx <- idx + 1L
      }
    }
  }
  fourier_conv <- bind_rows(fourier_rows) %>%
    group_by(tau_days, k) %>%
    mutate(
      reference_u_max = max(u_max),
      reference_iv = cf_lewis_iv[u_max == max(u_max)][1],
      reference_price = cf_lewis_price[u_max == max(u_max)][1],
      iv_error_to_largest_u = cf_lewis_iv - reference_iv,
      price_error_to_largest_u = cf_lewis_price - reference_price
    ) %>%
    ungroup()
  write_csv(fourier_conv, file.path(out_dir, "phaseS_rheston_fourier_truncation_convergence.csv"))
  print(fourier_conv)
}

# ============================================================
# 13. H-sensitivity of ATM skew
# ============================================================

if (RUN_H_SENSITIVITY) {
  cat("Running H-sensitivity check for ATM skew term structure...\n")
  sensitivity_k <- intersect(c(-0.05, 0, 0.05), unique(pricing_grid$k))
  if (length(sensitivity_k) < 3L) sensitivity_k <- unique(pricing_grid$k)[abs(unique(pricing_grid$k)) <= CORE_K_MAX + 1e-12]
  sensitivity_grid <- pricing_grid %>%
    filter(k %in% sensitivity_k) %>%
    arrange(tau_days, k)
  sens_surfaces <- list()
  idx <- 1L
  for (Hval in SENSITIVITY_H_GRID) {
    par_h <- rh_par
    par_h$alpha <- Hval + 0.5
    par_h$H <- Hval
    cat("  H sensitivity:", Hval, "alpha =", par_h$alpha, "\n")
    surf_h <- price_surface_from_par(par_h, sensitivity_grid, progress = FALSE)$prices %>%
      mutate(H = Hval, alpha = Hval + 0.5)
    sens_surfaces[[idx]] <- surf_h
    idx <- idx + 1L
  }
  h_sensitivity_surface <- bind_rows(sens_surfaces)
  h_sensitivity_skew <- h_sensitivity_surface %>%
    group_by(H, alpha, tau_days, tau) %>%
    group_modify(~ estimate_atm_skew(.x, "cf_lewis_iv", window = 0.055)) %>%
    ungroup() %>%
    rename(cf_atm_skew = psi, cf_atm_iv = sigma_atm, n_local = n_local)
  h_sensitivity_summary <- h_sensitivity_skew %>%
    filter(is.finite(cf_atm_skew), abs(cf_atm_skew) > 0) %>%
    group_by(H, alpha) %>%
    summarise(
      estimated_loglog_slope = coef(lm(log(abs(cf_atm_skew)) ~ log(tau)))[2],
      theoretical_H_minus_half = first(H) - 0.5,
      recovered_H_from_slope = estimated_loglog_slope + 0.5,
      n_maturities = n(),
      .groups = "drop"
    )
  write_csv(h_sensitivity_surface, file.path(out_dir, "phaseS_rheston_H_sensitivity_surface.csv"))
  write_csv(h_sensitivity_skew, file.path(out_dir, "phaseS_rheston_H_sensitivity_atm_skew.csv"))
  write_csv(h_sensitivity_summary, file.path(out_dir, "phaseS_rheston_H_sensitivity_summary.csv"))
  print(h_sensitivity_summary)
}

# ============================================================
# 14. Report-ready plots
# ============================================================

maturity_palette <- c("7d"="#1B9E77","14d"="#D95F02","30d"="#7570B3","90d"="#E7298A","180d"="#1F78B4")
method_palette <- c("rough Heston CF + Lewis"="#0072B2","iVi reference MC"="#D55E00")
region_palette <- c("Core: |k| <= 0.05"="#0072B2","Outer wing"="#999999","7d stress: k = 0.10"="#D55E00","Full grid"="#666666")

region_label_fun <- function(x) {
  dplyr::recode(
    x,
    headline_core_abs_k_le_0p05 = "Core: |k| <= 0.05",
    full_grid_all_points = "Full grid",
    stress_7d_right_wing_k_0p10 = "7d stress: k = 0.10",
    outer_wing_nonheadline = "Outer wing"
  )
}

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5, colour = "grey20"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
    plot.margin = margin(8, 10, 8, 10)
  )

cf_vs_ivi_plot <- cf_vs_ivi %>%
  mutate(
    tau_label = factor(paste0(tau_days, "d"), levels = paste0(sort(unique(tau_days)), "d")),
    validation_label = factor(region_label_fun(validation_region), levels = c("Core: |k| <= 0.05", "Outer wing", "7d stress: k = 0.10"))
  )

plot_cf <- rheston_cf_prices %>%
  transmute(tau_days, tau, k, iv = cf_lewis_iv, series = "rough Heston CF + Lewis", iv_se = NA_real_)
plot_ivi <- ivi_reference_surface %>%
  transmute(tau_days, tau, k, iv = mc_iv, series = "iVi reference MC", iv_se = mc_iv_se)
plot_iv_compare <- bind_rows(plot_cf, plot_ivi) %>%
  mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(sort(unique(tau_days)), "d")),
         series = factor(series, levels = c("rough Heston CF + Lewis", "iVi reference MC"))) %>%
  arrange(tau_days, series, k)

p_iv_compare <- ggplot(plot_iv_compare, aes(x = k, y = iv, colour = series, linetype = series, shape = series)) +
  geom_errorbar(data = plot_iv_compare %>% filter(series == "iVi reference MC"), aes(ymin = iv - 2 * iv_se, ymax = iv + 2 * iv_se), width = 0.004, linewidth = 0.25, alpha = 0.65) +
  geom_line(linewidth = 0.70, na.rm = TRUE) +
  geom_point(size = 1.55, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey35") +
  facet_wrap(~ tau_label, scales = "free_y") +
  scale_colour_manual(values = method_palette) +
  scale_linetype_manual(values = c("rough Heston CF + Lewis" = "solid", "iVi reference MC" = "dashed")) +
  scale_shape_manual(values = c("rough Heston CF + Lewis" = 16, "iVi reference MC" = 17)) +
  labs(title = "rough Heston CF/Lewis implied volatility validation", x = "Forward log-moneyness, k = log(K/F)", y = "B&S implied volatility", colour = NULL, linetype = NULL, shape = NULL) +
  plot_theme
plot_save(p_iv_compare, "plot_phaseS_rheston_cf_vs_ivi_iv_slices", width = 11, height = 6.2)

p_iv_errors <- cf_vs_ivi_plot %>%
  ggplot(aes(x = k, y = iv_error_cf_minus_ivi)) +
  geom_ribbon(aes(ymin = -2 * ivi_ref_iv_se, ymax = 2 * ivi_ref_iv_se), fill = "grey75", alpha = 0.45, colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey25") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey35") +
  geom_line(linewidth = 0.60, colour = "grey15", na.rm = TRUE) +
  geom_point(aes(colour = validation_label), size = 1.8, na.rm = TRUE) +
  facet_wrap(~ tau_label, scales = "free_y") +
  scale_colour_manual(values = region_palette) +
  labs(title = "IV error of rough Heston CF/Lewis against iVi reference", x = "Forward log-moneyness, k = log(K/F)", y = "CF IV - iVi reference IV", colour = "Reporting region") +
  plot_theme
plot_save(p_iv_errors, "plot_phaseS_rheston_cf_minus_ivi_iv_errors", width = 11, height = 6.2)

atm_plot_tbl <- atm_skew_comparison %>%
  select(tau_days, tau, `CF/Lewis` = cf_atm_skew, `iVi reference` = ivi_atm_skew) %>%
  pivot_longer(cols = c(`CF/Lewis`, `iVi reference`), names_to = "series", values_to = "atm_skew") %>%
  mutate(series = factor(series, levels = c("CF/Lewis", "iVi reference")))

p_atm_skew <- atm_plot_tbl %>%
  ggplot(aes(x = tau_days, y = atm_skew, colour = series, linetype = series, shape = series)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_line(linewidth = 0.75, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  labs(title = "ATM-skew comparison: CF/Lewis versus iVi", x = "Maturity in calendar days", y = expression(paste("ATM skew, ", psi(tau))), colour = NULL, linetype = NULL, shape = NULL) +
  plot_theme
plot_save(p_atm_skew, "plot_phaseS_rheston_cf_vs_ivi_atm_skew", width = 8.5, height = 5.2)

region_plot_tbl <- summary_global_report %>%
  mutate(validation_label = factor(region_label_fun(validation_region), levels = c("Core: |k| <= 0.05", "Full grid", "Outer wing", "7d stress: k = 0.10")))
region_ymax <- max(2.2, max(region_plot_tbl$iv_rmse_over_mean_mcse, na.rm = TRUE) * 1.20)

p_region_rmse <- region_plot_tbl %>%
  ggplot(aes(x = validation_label, y = iv_rmse_over_mean_mcse, fill = validation_label)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey35") +
  geom_hline(yintercept = 2, linewidth = 0.35, linetype = "dashed", colour = "grey35") +
  geom_text(aes(label = sprintf("%.2f", iv_rmse_over_mean_mcse)), hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = region_palette) +
  scale_y_continuous(limits = c(0, region_ymax), expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(title = "CF/Lewis validation metric by reporting region", x = NULL, y = "IV RMSE / mean iVi IV SE", fill = NULL) +
  plot_theme +
  theme(legend.position = "none")
plot_save(p_region_rmse, "plot_phaseS_rheston_cf_region_rmse_over_mcse", width = 9, height = 5.2)

if (RUN_H_SENSITIVITY && exists("h_sensitivity_summary")) {
  p_h_sens <- h_sensitivity_summary %>%
    ggplot(aes(x = H, y = estimated_loglog_slope)) +
    geom_abline(intercept = -0.5, slope = 1, linewidth = 0.55, linetype = "dashed", colour = "grey35") +
    geom_point(size = 2.2, colour = "#0072B2") +
    geom_line(linewidth = 0.75, colour = "#0072B2") +
    labs(title = "H-sensitivity of the ATM-skew slope", x = "Hurst parameter H", y = expression(paste("Estimated slope of log |", psi(tau), "|"))) +
    plot_theme
  plot_save(p_h_sens, "plot_phaseS_rheston_H_sensitivity_atm_skew_slope", width = 8.5, height = 5.2)
}

# ============================================================
# 15. Figure manifest and final summary
# ============================================================

figure_files <- c(
  "plot_phaseS_rheston_cf_vs_ivi_iv_slices",
  "plot_phaseS_rheston_cf_minus_ivi_iv_errors",
  "plot_phaseS_rheston_cf_vs_ivi_atm_skew",
  "plot_phaseS_rheston_cf_region_rmse_over_mcse"
)
figure_use <- c(TRUE, TRUE, TRUE, TRUE)
figure_desc <- c(
  "rough Heston CF/Lewis implied-volatility slices compared with iVi Monte Carlo reference.",
  "Implied volatility error CF minus iVi reference with Monte Carlo error band and reporting-region labels.",
  "ATM-skew comparison between rough Heston CF/Lewis and iVi reference.",
  "Global validation metrics separated by headline core region, full grid, stress point, and outer wing."
)
if (RUN_H_SENSITIVITY) {
  figure_files <- c(figure_files, "plot_phaseS_rheston_H_sensitivity_atm_skew_slope")
  figure_use <- c(figure_use, TRUE)
  figure_desc <- c(figure_desc, "Sensitivity of the estimated ATM-skew log-log slope to the Hurst parameter.")
}

figure_manifest <- tibble(file = paste0(figure_files, ".", PLOT_FORMAT), use_in_thesis = figure_use, description = figure_desc)
write_csv(figure_manifest, file.path(out_dir, "phaseS_rheston_cf_pricing_figure_manifest.csv"))

bad_cf_count <- sum(!cf_eval_grid$cf_ok | !is.finite(cf_eval_grid$cf_real) | !is.finite(cf_eval_grid$cf_imag))
price_clip_count <- sum(rheston_cf_prices$cf_price_clipped, na.rm = TRUE)
na_iv_count <- sum(!is.finite(rheston_cf_prices$cf_lewis_iv))

headline <- summary_global_core
full <- summary_global_full
stress <- summary_global_stress

# Final pass/fail table for pricing-engine validation.
acceptance_checks <- tibble(
  check = c(
    "B&S/Lewis max abs price error < 1e-4",
    "max |phi(0)-1| < 1e-8",
    "max |phi(-i)-1| < 1e-8",
    "headline core clipping count = 0",
    "headline IV RMSE <= 2 * mean iVi IV SE",
    "ATM skew extraction uses two-sided local-quadratic support"
  ),
  value = c(
    black_sanity$summary$max_abs_error[1],
    max(martingale_checks$abs_cf0_minus_1, na.rm = TRUE),
    max(martingale_checks$abs_cf_minus_i_minus_1, na.rm = TRUE),
    core_clip_count,
    headline$iv_rmse_over_mean_mcse[1],
    sum(atm_skew_extraction_audit$headline_skew_estimator_ok, na.rm = TRUE)
  ),
  pass = c(
    black_sanity$summary$max_abs_error[1] < 1e-4,
    max(martingale_checks$abs_cf0_minus_1, na.rm = TRUE) < 1e-8,
    max(martingale_checks$abs_cf_minus_i_minus_1, na.rm = TRUE) < 1e-8,
    core_clip_count == 0,
    headline$iv_rmse_over_mean_mcse[1] <= 2,
    isTRUE(all(atm_skew_extraction_audit$headline_skew_estimator_ok))
  )
)
write_csv(acceptance_checks, file.path(out_dir, "phaseS_rheston_cf_acceptance_checks.csv"))
print(acceptance_checks)

cat("\n=== Phase S rough-Heston CF pricing completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("FAST_TEST mode:", FAST_TEST, "\n")
cat("Maturities priced:", paste(maturity_list, collapse = ", "), "days\n")
cat("Number of strikes/log-moneyness points:", length(unique(pricing_grid$k)), "\n")
cat("DU target:", DU_TARGET, "\n")
cat("Riccati corrector iterations:", RICCATI_CORRECTOR_ITERS, "\n")
cat("Bad CF node count:", bad_cf_count, "\n")
cat("CF price clipping count:", price_clip_count, "\n")
cat("Headline core clipping count:", core_clip_count, "\n")
cat("CF NA IV count:", na_iv_count, "\n\n")

cat("=== Headline ATM/near-ATM validation: |k| <= ", CORE_K_MAX, " ===\n", sep = "")
cat("Headline quotes:", headline$n_quotes[1], "\n")
cat("Headline IV RMSE:", signif(headline$iv_rmse[1], 5), "\n")
cat("Headline IV MAE:", signif(headline$iv_mae[1], 5), "\n")
cat("Headline IV RMSE / mean iVi MCSE:", signif(headline$iv_rmse_over_mean_mcse[1], 5), "\n")
cat("Headline price RMSE:", signif(headline$price_rmse[1], 5), "\n")
cat("Headline price RMSE / mean iVi price MCSE:", signif(headline$price_rmse_over_mean_mcse[1], 5), "\n\n")

cat("=== Full-grid diagnostic, not headline ===\n")
cat("Full-grid quotes:", full$n_quotes[1], "\n")
cat("Full-grid IV RMSE:", signif(full$iv_rmse[1], 5), "\n")
cat("Full-grid IV RMSE / mean iVi MCSE:", signif(full$iv_rmse_over_mean_mcse[1], 5), "\n\n")

if (nrow(stress) > 0 && is.finite(stress$n_quotes[1]) && stress$n_quotes[1] > 0) {
  cat("=== Stress point: 7d, k = 0.10 ===\n")
  cat("Stress quotes:", stress$n_quotes[1], "\n")
  cat("Stress IV RMSE:", signif(stress$iv_rmse[1], 5), "\n")
  cat("Stress IV RMSE / mean iVi MCSE:", signif(stress$iv_rmse_over_mean_mcse[1], 5), "\n")
  cat("Stress max abs price error:", signif(stress$max_abs_price_error[1], 5), "\n")
  cat("Stress max abs price error / MCSE:", signif(stress$max_abs_price_error_over_mcse[1], 5), "\n\n")
}

cat("Main files to inspect:\n")
cat("  phaseS_rheston_cf_acceptance_checks.csv\n")
cat("  phaseS_rheston_cf_vs_ivi_summary_global_report.csv\n")
cat("  phaseS_rheston_cf_vs_ivi_atm_skew_comparison.csv\n")
cat("  phaseS_rheston_atm_skew_extraction_audit.csv\n")
cat("  phaseS_rheston_atm_skew_central_difference_check.csv\n")
cat("  phaseS_rheston_atm_skew_loglog_summary.csv\n")
cat("  phaseS_alpha_one_heston_limit_check.csv\n")
cat("  phaseS_rheston_cf_grid_summary.csv\n")
cat("Files written: CSV tables and PDF plots.\n")
