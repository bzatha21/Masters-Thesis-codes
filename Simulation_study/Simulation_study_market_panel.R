#!/usr/bin/env Rscript

# ============================================================
# Phase S: controlled simulation study
# ------------------------------------------------------------
#
# Purpose:
#   Construct two controlled validation layers before empirical
#   calibration:
#
#   Experiment 1: synthetic SPX-scale IV surface with imposed rough-skew
#                 scaling.
#   Experiment 2: iVi Volterra Heston Monte Carlo reference surface.
#
# Scope discipline:
#   - Experiment 1 checks the ATM-skew extraction and recovery workflow.
#   - The H = 1/2 comparator is a smooth reduced-form benchmark.
#   - Experiment 2 generates an independent model-based reference surface.
#   - Rough Heston CF/Lewis pricing is handled separately in
#     Phase_s_rheston_cf_pricing.R.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(271828)

# ============================================================
# 0. Run configuration and output paths
# ============================================================

FAST_TEST <- identical(Sys.getenv("PHASES_FAST_TEST"), "1")

CLEAN_OUTPUT_DIR <- TRUE
PLOT_FORMAT <- "pdf"
SHOW_PLOTS <- !FAST_TEST
SAVE_PLOTS <- TRUE
RUN_IVI_MODEL_EXPERIMENT <- TRUE

# Monte Carlo grid settings; automatically reduced in FAST_TEST mode.
if (FAST_TEST) {
  IVI_N_PATHS <- 2000L
  IVI_REF_N_PATHS <- 4000L
  IVI_N_STEPS_GRID <- c(4L, 8L)
  IVI_REF_N_STEPS <- 12L
  IVI_TAU_DAYS <- c(7L, 30L)
  IVI_K_GRID <- c(-0.05, 0.00, 0.05)
} else {
  IVI_N_PATHS <- 100000L
  IVI_REF_N_PATHS <- 200000L
  IVI_N_STEPS_GRID <- c(8L, 16L, 32L, 64L)
  IVI_REF_N_STEPS <- 96L
  IVI_TAU_DAYS <- c(7L, 14L, 30L, 90L, 180L)
  IVI_K_GRID <- c(-0.10, -0.05, -0.025, 0.00, 0.025, 0.05, 0.10)
}

# SPX-scale anchor values used for the synthetic forward and strike grids.
S0 <- 1227.73
F0 <- S0
D0 <- 1
valuation_date <- as.Date("2005-09-15")
q_div <- 0.018

rate_curve <- tibble(
  tau_days = c(7, 14, 30, 60, 90, 180, 365),
  r = c(0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041)
)

# Interpolated deterministic rate used only to build synthetic forwards.
rate_for_days <- function(days) {
  approx(rate_curve$tau_days, rate_curve$r, xout = days, rule = 2)$y
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

out_dir <- file.path(base_dir, "data", "SimulationStudy", "phaseS_v2_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_OUTPUT_DIR) {
  old_files <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)
}

cat("Output directory:\n", out_dir, "\n\n", sep = "")
cat("FAST_TEST mode:", FAST_TEST, "\n\n")

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

# Persist the synthetic market conventions used by all downstream tables.
market_assumptions <- tibble(
  valuation_date = valuation_date,
  S0 = S0,
  q_div = q_div,
  F0 = F0,
  D0 = D0,
  fast_test = FAST_TEST,
  note = "Synthetic SPX-scale assumptions; not an empirical SPX option-chain calibration"
)

write_csv(market_assumptions, file.path(out_dir, "simulation_v2_market_assumptions.csv"))
write_csv(rate_curve, file.path(out_dir, "simulation_v2_synthetic_rate_curve.csv"))

# ============================================================
# 1. Parameter scenarios and interpretation notes
# ============================================================

# Published-calibration scale anchor; the iVi base scenario is specified separately.
parameter_scenarios <- tibble::tribble(
  ~scenario_id, ~scenario_type, ~gamma, ~theta, ~nu, ~rho, ~V0, ~alpha, ~use_in_main,
  "ER2019_SPX_anchor", "published_calibration_anchor", 0.100, 0.3156, 0.331, -0.681, 0.0392, 0.62, TRUE,
  "iVi_stable_base", "model_based_simulation_benchmark", 1.500, 0.0400, 0.450, -0.700, 0.0450, 0.60, TRUE,
  "H_0p05_sensitivity", "roughness_sensitivity", 1.500, 0.0400, 0.450, -0.700, 0.0450, 0.55, FALSE,
  "H_0p20_sensitivity", "roughness_sensitivity", 1.500, 0.0400, 0.450, -0.700, 0.0450, 0.70, FALSE,
  "H_0p50_smooth_limit", "smoothness_sensitivity", 1.500, 0.0400, 0.450, -0.700, 0.0450, 1.00, FALSE,
  "rho_low_sensitivity", "skew_sensitivity", 1.500, 0.0400, 0.450, -0.500, 0.0450, 0.60, FALSE,
  "rho_high_sensitivity", "skew_sensitivity", 1.500, 0.0400, 0.450, -0.850, 0.0450, 0.60, FALSE,
  "nu_low_sensitivity", "smile_sensitivity", 1.500, 0.0400, 0.300, -0.700, 0.0450, 0.60, FALSE,
  "nu_high_sensitivity", "smile_sensitivity", 1.500, 0.0400, 0.600, -0.700, 0.0450, 0.60, FALSE
) %>%
  mutate(H = alpha - 0.5)

write_csv(parameter_scenarios, file.path(out_dir, "simulation_v2_parameter_scenarios.csv"))
print(parameter_scenarios)

# ============================================================
# 2. Black-Scholes functions in forward form
# ============================================================

# Prices are expressed as D_tau times a forward-price expectation.
black_call_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) return(D * max(F - K, 0))
  vol_sqrt <- sigma * sqrt(tau)
  d1 <- (log(F / K) + 0.5 * sigma^2 * tau) / vol_sqrt
  d2 <- d1 - vol_sqrt
  D * (F * pnorm(d1) - K * pnorm(d2))
}

black_put_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) return(D * max(K - F, 0))
  vol_sqrt <- sigma * sqrt(tau)
  d1 <- (log(F / K) + 0.5 * sigma^2 * tau) / vol_sqrt
  d2 <- d1 - vol_sqrt
  D * (K * pnorm(-d2) - F * pnorm(-d1))
}

black_price_forward <- function(cp, F, K, tau, sigma, D = 1) {
  if (cp == "C") black_call_forward(F, K, tau, sigma, D)
  else if (cp == "P") black_put_forward(F, K, tau, sigma, D)
  else stop("cp must be 'C' or 'P'.")
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
# 3. Black-Scholes/Lewis control check
# ============================================================

# Characteristic function of the log-forward return under Black--Scholes.
black_cf_log_forward <- function(z, tau, sigma) {
  exp(-0.5 * sigma^2 * tau * (z^2 + 1i * z))
}

lewis_call_from_cf <- function(k, tau, cf_fun, F = F0, D = D0, u_max = 100, n_u = 2000) {
  K <- F * exp(k)
  u <- seq(0, u_max, length.out = n_u + 1L)
  du <- u_max / n_u
  w <- rep(du, length(u))
  w[1L] <- 0.5 * du
  w[length(w)] <- 0.5 * du
  cf_vec <- sapply(u, function(x) cf_fun(x - 0.5i))
  vals <- exp(-1i * u * k) * cf_vec / (u^2 + 0.25)
  D * (F - sqrt(F * K) / pi * sum(w * Re(vals)))
}

run_black_lewis_sanity <- function() {
  test_grid <- expand.grid(k = c(-0.15, 0, 0.15), tau = c(30, 180) / 365)
  sigma <- 0.22
  err <- numeric(nrow(test_grid))
  for (i in seq_len(nrow(test_grid))) {
    k_i <- test_grid$k[i]
    tau_i <- test_grid$tau[i]
    K_i <- F0 * exp(k_i)
    p_black <- black_call_forward(F0, K_i, tau_i, sigma, D0)
    p_lewis <- lewis_call_from_cf(
      k = k_i,
      tau = tau_i,
      cf_fun = function(z) black_cf_log_forward(z, tau_i, sigma),
      F = F0,
      D = D0,
      u_max = 100,
      n_u = if (FAST_TEST) 800 else 2000
    )
    err[i] <- abs(p_black - p_lewis)
  }
  tibble(max_abs_error = max(err), pass = max(err) < 1e-4)
}

black_test <- run_black_lewis_sanity()
write_csv(black_test, file.path(out_dir, "simulation_v2_black_lewis_sanity.csv"))
print(black_test)
if (!isTRUE(black_test$pass[1])) stop("Black/Lewis sanity check failed.")

# ============================================================
# 4. Experiment 1: controlled rough-skew target surface
# ============================================================

rough_design <- list(
  H = 0.10,
  alpha = 0.60,
  sigma_long = 0.180,
  sigma_short_amp = 0.040,
  sigma_decay = 0.35,
  A_skew = 0.035,
  curvature0 = 0.006,
  curvature1 = 0.008,
  eta = 0.055
)

sigma_atm_fun <- function(tau, p = rough_design) {
  p$sigma_long + p$sigma_short_amp * exp(-tau / p$sigma_decay)
}

psi_rough_fun <- function(tau, p = rough_design) {
  -p$A_skew * tau^(p$H - 0.5)
}

curvature_fun <- function(tau, p = rough_design) {
  p$curvature0 + p$curvature1 * sqrt(tau)
}

# Synthetic total-variance surface with imposed rough-skew power law.
rough_surface_iv <- function(k, tau, p = rough_design) {
  sig0 <- sigma_atm_fun(tau, p)
  psi <- psi_rough_fun(tau, p)
  w0 <- sig0^2 * tau
  w_slope <- 2 * tau * sig0 * psi
  curv <- curvature_fun(tau, p)
  w <- w0 + w_slope * k + curv * (sqrt(k^2 + p$eta^2) - p$eta)
  sqrt(pmax(w, 1e-10) / tau)
}

rough_design_tbl <- tibble(
  object = "controlled_rough_skew_synthetic_surface",
  H = rough_design$H,
  alpha = rough_design$alpha,
  sigma_long = rough_design$sigma_long,
  sigma_short_amp = rough_design$sigma_short_amp,
  sigma_decay = rough_design$sigma_decay,
  A_skew = rough_design$A_skew,
  curvature0 = rough_design$curvature0,
  curvature1 = rough_design$curvature1,
  eta = rough_design$eta,
  reference_slope = rough_design$H - 0.5
)
write_csv(rough_design_tbl, file.path(out_dir, "simulation_v2_design_parameters.csv"))
print(rough_design_tbl)

# Build maturity-specific forward log-moneyness grids, k = log(K/F_tau).
make_market_grid <- function(tau_days_vec, S0, q_div, strike_step = 5) {
  bind_rows(lapply(tau_days_vec, function(td) {
    tau <- td / 365
    r_tau <- rate_for_days(td)
    D_tau <- exp(-r_tau * tau)
    F_tau <- S0 * exp((r_tau - q_div) * tau)
    K_min <- floor(F_tau * exp(-0.16) / strike_step) * strike_step
    K_max <- ceiling(F_tau * exp(0.12) / strike_step) * strike_step
    tibble(
      tau_days = td,
      tau = tau,
      r = r_tau,
      q = q_div,
      D = D_tau,
      F = F_tau,
      K = seq(K_min, K_max, by = strike_step)
    ) %>% mutate(k = log(K / F))
  }))
}

tau_days_vec <- if (FAST_TEST) c(7, 30, 90, 180) else c(7, 14, 21, 30, 45, 60, 90, 120, 180, 270, 365)
market_grid <- make_market_grid(tau_days_vec, S0, q_div, strike_step = 5)

maturity_shocks <- tibble(
  tau_days = tau_days_vec,
  mat_noise = rnorm(length(tau_days_vec), mean = 0, sd = 0.0015),
  skew_noise = rnorm(length(tau_days_vec), mean = 0, sd = 0.0040),
  curv_noise = rnorm(length(tau_days_vec), mean = 0, sd = 0.0030)
)

# Add quote-like perturbations while preserving the controlled rough-skew design.
surface_grid <- market_grid %>%
  left_join(maturity_shocks, by = "tau_days") %>%
  group_by(tau_days) %>%
  mutate(k2_centered = k^2 - mean(k^2)) %>%
  ungroup() %>%
  mutate(
    target_iv_clean = mapply(rough_surface_iv, k, tau),
    residual_iv = mat_noise + skew_noise * k + curv_noise * k2_centered +
      0.0035 * exp(-tau / 0.15) * exp(-0.5 * ((k + 0.115) / 0.025)^2) -
      0.0015 * exp(-tau / 0.30) * exp(-0.5 * ((k - 0.075) / 0.035)^2),
    target_iv = pmax(target_iv_clean + residual_iv, 0.03),
    target_call_price = mapply(function(F_, K_, tau_, iv_, D_) black_call_forward(F_, K_, tau_, iv_, D_), F, K, tau, target_iv, D),
    target_put_price = mapply(function(F_, K_, tau_, iv_, D_) black_put_forward(F_, K_, tau_, iv_, D_), F, K, tau, target_iv, D),
    cp = ifelse(K >= F, "C", "P"),
    target_price = ifelse(cp == "C", target_call_price, target_put_price),
    intrinsic = ifelse(cp == "C", D * pmax(F - K, 0), D * pmax(K - F, 0)),
    half_spread_iv = 0.0015 + 0.0060 * abs(k) + 0.0030 / sqrt(pmax(tau_days, 7) / 30),
    micro_noise_sd = 0.20 * half_spread_iv,
    obs_iv_mid = pmax(target_iv + rnorm(n(), mean = 0, sd = micro_noise_sd), 0.03),
    bid_iv = pmax(obs_iv_mid - half_spread_iv, 0.01),
    ask_iv = pmax(obs_iv_mid + half_spread_iv, bid_iv + 0.0005),
    obs_iv = 0.5 * (bid_iv + ask_iv),
    bid_price = mapply(function(cp_, F_, K_, tau_, iv_, D_) black_price_forward(cp_, F_, K_, tau_, iv_, D_), cp, F, K, tau, bid_iv, D),
    ask_price = mapply(function(cp_, F_, K_, tau_, iv_, D_) black_price_forward(cp_, F_, K_, tau_, iv_, D_), cp, F, K, tau, ask_iv, D),
    obs_price = 0.5 * (bid_price + ask_price),
    spread_iv = ask_iv - bid_iv,
    weight = 1 / pmax(spread_iv^2, 1e-7),
    valid_target = is.finite(target_iv) & is.finite(target_price) & target_iv > 0 & target_price >= intrinsic,
    tau_label = paste0(tau_days, "d")
  )

if (any(!surface_grid$valid_target)) {
  write_csv(surface_grid %>% filter(!valid_target), file.path(out_dir, "simulation_v2_invalid_target_rows.csv"))
  stop("Synthetic target surface contains invalid rows.")
}

surface_static_checks <- surface_grid %>%
  mutate(total_variance = target_iv^2 * tau) %>%
  group_by(tau_days, tau) %>%
  summarise(
    min_total_variance = min(total_variance),
    min_iv = min(target_iv),
    max_iv = max(target_iv),
    min_price_minus_intrinsic = min(target_price - intrinsic),
    .groups = "drop"
  )

write_csv(surface_grid, file.path(out_dir, "simulation_v2_synthetic_surface.csv"))
write_csv(surface_static_checks, file.path(out_dir, "simulation_v2_surface_static_checks.csv"))
cat("Synthetic surface rows:", nrow(surface_grid), "\n")

# ============================================================
# 5. ATM skew extraction and reduced-form recovery
# ============================================================

estimate_atm_skew <- function(df, iv_col = "target_iv", window = 0.055) {
  # ATM skew is the local derivative of the IV slice at k = 0.
  # The controlled grid requires a local quadratic fit with near-ATM support.
  d <- df %>% filter(abs(k) <= window + 1e-12, is.finite(.data[[iv_col]]))
  if (nrow(d) < 5L) {
    return(tibble(
      psi = NA_real_,
      intercept = NA_real_,
      n_local = nrow(d),
      fit_type = "insufficient"
    ))
  }
  d <- d %>% mutate(local_weight = pmax(0, 1 - abs(k) / window))
  fit <- lm(as.formula(paste0(iv_col, " ~ k + I(k^2)")), data = d, weights = local_weight)
  cf <- coef(fit)
  tibble(
    psi = unname(cf["k"]),
    intercept = unname(cf["(Intercept)"]),
    n_local = nrow(d),
    fit_type = "local_quadratic"
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

atm_skew_tbl <- surface_grid %>%
  group_by(tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "target_iv", window = 0.055)) %>%
  ungroup() %>%
  rename(
    psi_target_est = psi,
    sigma_atm_target_est = intercept,
    n_local_target = n_local,
    fit_type_target = fit_type
  ) %>%
  left_join(
    surface_grid %>%
      group_by(tau_days, tau) %>%
      group_modify(~ estimate_atm_skew(.x, "obs_iv", window = 0.055)) %>%
      ungroup() %>%
      rename(
        psi_obs_est = psi,
        sigma_atm_obs_est = intercept,
        n_local_obs = n_local,
        fit_type_obs = fit_type
      ),
    by = c("tau_days", "tau")
  ) %>%
  mutate(analytic_psi = psi_rough_fun(tau), abs_psi = abs(psi_target_est), reference_power = rough_design$H - 0.5)

slope_fit <- lm(log(abs_psi) ~ log(tau), data = atm_skew_tbl %>% filter(is.finite(abs_psi), abs_psi > 0))
atm_slope <- unname(coef(slope_fit)[2])
atm_intercept <- unname(coef(slope_fit)[1])
atm_summary <- tibble(
  estimated_loglog_slope = atm_slope,
  reference_H_minus_half = rough_design$H - 0.5,
  H_recovered_from_slope = atm_slope + 0.5,
  n_maturities = nrow(atm_skew_tbl),
  passes_H_tolerance_0p03 = abs((atm_slope + 0.5) - rough_design$H) < 0.03
)

synthetic_near_atm_support <- surface_grid %>%
  group_by(tau_days, tau) %>%
  group_modify(~ near_atm_support_audit(.x, window = 0.055)) %>%
  ungroup()

synthetic_cd_target_h0025 <- surface_grid %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "target_iv", h = 0.025)) %>%
  ungroup() %>%
  rename(target_cd_skew_h0025 = cd_skew, target_cd_available_h0025 = cd_available)

synthetic_cd_obs_h0025 <- surface_grid %>%
  group_by(tau_days, tau) %>%
  group_modify(~ central_difference_skew(.x, "obs_iv", h = 0.025)) %>%
  ungroup() %>%
  rename(obs_cd_skew_h0025 = cd_skew, obs_cd_available_h0025 = cd_available)

atm_skew_extraction_audit <- atm_skew_tbl %>%
  left_join(synthetic_near_atm_support, by = c("tau_days", "tau")) %>%
  left_join(synthetic_cd_target_h0025, by = c("tau_days", "tau")) %>%
  left_join(synthetic_cd_obs_h0025, by = c("tau_days", "tau")) %>%
  mutate(
    headline_estimator_target_ok = fit_type_target == "local_quadratic" & n_local_target >= 5L,
    headline_estimator_obs_ok = fit_type_obs == "local_quadratic" & n_local_obs >= 5L,
    headline_skew_estimator_ok = two_sided_support & headline_estimator_target_ok & headline_estimator_obs_ok,
    target_lq_minus_cd_h0025 = psi_target_est - target_cd_skew_h0025,
    obs_lq_minus_cd_h0025 = psi_obs_est - obs_cd_skew_h0025
  )

atm_loglog_points <- atm_skew_tbl %>%
  filter(is.finite(abs_psi), abs_psi > 0) %>%
  mutate(
    log_abs_psi = log(abs_psi),
    fitted_log_abs_psi = fitted(slope_fit),
    loglog_residual = residuals(slope_fit)
  ) %>%
  select(tau_days, tau, psi_target_est, abs_psi, log_abs_psi, fitted_log_abs_psi, loglog_residual)

write_csv(atm_skew_tbl, file.path(out_dir, "simulation_v2_atm_skew.csv"))
write_csv(atm_summary, file.path(out_dir, "simulation_v2_atm_skew_summary.csv"))
write_csv(atm_skew_extraction_audit, file.path(out_dir, "simulation_v2_atm_skew_extraction_audit.csv"))
write_csv(atm_loglog_points, file.path(out_dir, "simulation_v2_atm_skew_loglog_points.csv"))
print(atm_summary)
print(atm_skew_extraction_audit)

# Reduced-form parameter recovery with a held-out split; 
logit <- function(x, lo, hi) log((x - lo) / (hi - x))
inv_logit <- function(y, lo, hi) lo + (hi - lo) / (1 + exp(-y))

rough_bounds <- list(
  sigma_long = c(0.12, 0.26), sigma_short_amp = c(0.00, 0.10), sigma_decay = c(0.05, 1.50),
  A_skew = c(0.005, 0.080), H = c(0.03, 0.30), curvature0 = c(0.000, 0.030), curvature1 = c(0.000, 0.040)
)

pack_rough <- function(p) c(
  logit(p$sigma_long, rough_bounds$sigma_long[1], rough_bounds$sigma_long[2]),
  logit(p$sigma_short_amp, rough_bounds$sigma_short_amp[1], rough_bounds$sigma_short_amp[2]),
  logit(p$sigma_decay, rough_bounds$sigma_decay[1], rough_bounds$sigma_decay[2]),
  logit(p$A_skew, rough_bounds$A_skew[1], rough_bounds$A_skew[2]),
  logit(p$H, rough_bounds$H[1], rough_bounds$H[2]),
  logit(p$curvature0, rough_bounds$curvature0[1], rough_bounds$curvature0[2]),
  logit(p$curvature1, rough_bounds$curvature1[1], rough_bounds$curvature1[2])
)

unpack_rough <- function(x) {
  H_value <- inv_logit(x[5], rough_bounds$H[1], rough_bounds$H[2])
  list(
    sigma_long = inv_logit(x[1], rough_bounds$sigma_long[1], rough_bounds$sigma_long[2]),
    sigma_short_amp = inv_logit(x[2], rough_bounds$sigma_short_amp[1], rough_bounds$sigma_short_amp[2]),
    sigma_decay = inv_logit(x[3], rough_bounds$sigma_decay[1], rough_bounds$sigma_decay[2]),
    A_skew = inv_logit(x[4], rough_bounds$A_skew[1], rough_bounds$A_skew[2]),
    H = H_value,
    alpha = H_value + 0.5,
    curvature0 = inv_logit(x[6], rough_bounds$curvature0[1], rough_bounds$curvature0[2]),
    curvature1 = inv_logit(x[7], rough_bounds$curvature1[1], rough_bounds$curvature1[2]),
    eta = rough_design$eta
  )
}

surface_iv_parametric <- function(k, tau, p, H_override = NULL) {
  H_use <- if (is.null(H_override)) p$H else H_override
  sig0 <- p$sigma_long + p$sigma_short_amp * exp(-tau / p$sigma_decay)
  psi <- -p$A_skew * tau^(H_use - 0.5)
  w0 <- sig0^2 * tau
  w_slope <- 2 * tau * sig0 * psi
  curv <- p$curvature0 + p$curvature1 * sqrt(tau)
  w <- w0 + w_slope * k + curv * (sqrt(k^2 + p$eta^2) - p$eta)
  sqrt(pmax(w, 1e-10) / tau)
}

set.seed(314159)
surface_grid <- surface_grid %>% group_by(tau_days) %>% arrange(k, .by_group = TRUE) %>% mutate(calib_flag = row_number() %% 3 != 0) %>% ungroup()
calib_panel <- surface_grid %>% filter(k >= -0.13, k <= 0.09, calib_flag)
test_panel <- surface_grid %>% filter(k >= -0.13, k <= 0.09, !calib_flag)

# Objective for the reduced-form rough-skew surface fit.
rough_obj <- function(x) {
  p <- unpack_rough(x)
  pred <- mapply(function(k_, tau_) surface_iv_parametric(k_, tau_, p), calib_panel$k, calib_panel$tau)
  if (any(!is.finite(pred))) return(1e8)
  mean(calib_panel$weight * (pred - calib_panel$obs_iv)^2)
}

# Same surface family with H fixed at the smooth benchmark value 1/2.
smooth_obj <- function(x) {
  p <- unpack_rough(x)
  pred <- mapply(function(k_, tau_) surface_iv_parametric(k_, tau_, p, H_override = 0.5), calib_panel$k, calib_panel$tau)
  if (any(!is.finite(pred))) return(1e8)
  mean(calib_panel$weight * (pred - calib_panel$obs_iv)^2)
}

x0 <- pack_rough(rough_design)
starts <- list(
  x0,
  x0 + c(0.20, -0.10, 0.15, 0.10, -0.25, 0.10, -0.10),
  x0 + c(-0.20, 0.15, -0.15, -0.10, 0.25, -0.10, 0.10),
  x0 + c(0.05, 0.10, -0.20, 0.20, -0.10, 0.05, 0.05),
  x0 + c(-0.10, -0.10, 0.25, 0.15, 0.15, 0.10, -0.20)
)

if (FAST_TEST) starts <- starts[1:2]

cat("Running reduced-form rough-skew recovery...\n")
rough_fits <- lapply(starts, function(s) optim(par = s, fn = rough_obj, method = "Nelder-Mead", control = list(maxit = if (FAST_TEST) 250 else 1600)))
rough_best <- rough_fits[[which.min(sapply(rough_fits, `[[`, "value"))]]
rough_hat <- unpack_rough(rough_best$par)

cat("Running smooth H=1/2 reduced-form benchmark...\n")
smooth_fits <- lapply(starts, function(s) optim(par = s, fn = smooth_obj, method = "Nelder-Mead", control = list(maxit = if (FAST_TEST) 250 else 1600)))
smooth_best <- smooth_fits[[which.min(sapply(smooth_fits, `[[`, "value"))]]
smooth_hat <- unpack_rough(smooth_best$par)

param_recovery <- bind_rows(
  tibble(model = "true_target_base_parameters", sigma_long = rough_design$sigma_long, sigma_short_amp = rough_design$sigma_short_amp, sigma_decay = rough_design$sigma_decay, A_skew = rough_design$A_skew, H = rough_design$H, alpha = rough_design$alpha, curvature0 = rough_design$curvature0, curvature1 = rough_design$curvature1, objective = NA_real_, convergence = NA_integer_),
  tibble(model = "rough_skew_recovered", sigma_long = rough_hat$sigma_long, sigma_short_amp = rough_hat$sigma_short_amp, sigma_decay = rough_hat$sigma_decay, A_skew = rough_hat$A_skew, H = rough_hat$H, alpha = rough_hat$alpha, curvature0 = rough_hat$curvature0, curvature1 = rough_hat$curvature1, objective = rough_best$value, convergence = rough_best$convergence),
  tibble(model = "smooth_H_0p5_reduced_form_benchmark", sigma_long = smooth_hat$sigma_long, sigma_short_amp = smooth_hat$sigma_short_amp, sigma_decay = smooth_hat$sigma_decay, A_skew = smooth_hat$A_skew, H = 0.5, alpha = 1.0, curvature0 = smooth_hat$curvature0, curvature1 = smooth_hat$curvature1, objective = smooth_best$value, convergence = smooth_best$convergence)
)

fit_panel <- surface_grid %>%
  mutate(
    rough_fit_iv = mapply(function(k_, tau_) surface_iv_parametric(k_, tau_, rough_hat), k, tau),
    smooth_fit_iv = mapply(function(k_, tau_) surface_iv_parametric(k_, tau_, smooth_hat, H_override = 0.5), k, tau),
    rough_fit_error_target = rough_fit_iv - target_iv,
    smooth_fit_error_target = smooth_fit_iv - target_iv,
    rough_fit_error_obs = rough_fit_iv - obs_iv,
    smooth_fit_error_obs = smooth_fit_iv - obs_iv
  )

error_summary <- fit_panel %>%
  filter(k >= -0.13, k <= 0.09, !calib_flag) %>%
  summarise(
    rough_iv_rmse_target = sqrt(mean(rough_fit_error_target^2)),
    smooth_iv_rmse_target = sqrt(mean(smooth_fit_error_target^2)),
    rough_iv_mae_target = mean(abs(rough_fit_error_target)),
    smooth_iv_mae_target = mean(abs(smooth_fit_error_target)),
    rough_iv_rmse_obs = sqrt(mean(rough_fit_error_obs^2)),
    smooth_iv_rmse_obs = sqrt(mean(smooth_fit_error_obs^2)),
    n_test_quotes = n(),
    .groups = "drop"
  ) %>%
  mutate(
    comparison = "held_out_rough_skew_recovery_vs_smooth_H_0p5_benchmark",
    rmse_ratio_smooth_over_rough_target = smooth_iv_rmse_target / rough_iv_rmse_target,
    rmse_ratio_smooth_over_rough_obs = smooth_iv_rmse_obs / rough_iv_rmse_obs
  )

write_csv(surface_grid, file.path(out_dir, "simulation_v2_synthetic_surface.csv"))
write_csv(param_recovery, file.path(out_dir, "simulation_v2_surface_model_parameter_recovery.csv"))
write_csv(fit_panel, file.path(out_dir, "simulation_v2_fit_panel_reduced_form.csv"))
write_csv(error_summary, file.path(out_dir, "simulation_v2_model_error_summary.csv"))
print(param_recovery)
print(error_summary)

# ============================================================
# 6. Simulated rough log-volatility diagnostic
# ============================================================

# Fractional Gaussian-noise path used only for the roughness diagnostic.
simulate_fgn_chol <- function(n, H) {
  h <- 0:(n - 1L)
  gamma_h <- 0.5 * (abs(h + 1)^(2 * H) - 2 * abs(h)^(2 * H) + abs(h - 1)^(2 * H))
  Sigma <- toeplitz(gamma_h)
  R <- chol(Sigma + diag(1e-10, n))
  as.numeric(t(R) %*% rnorm(n))
}

n_days_rough <- if (FAST_TEST) 350L else 1500L
dt_rough <- 1 / 252
H_true_rough <- rough_design$H
fgn <- simulate_fgn_chol(n_days_rough, H_true_rough) * dt_rough^H_true_rough
fbm_path <- cumsum(fgn)
rough_logvol_scale <- 0.15
logvol_path <- log(0.20) + rough_logvol_scale * fbm_path
vol_path <- exp(logvol_path)

rough_path_tbl <- tibble(day = seq_len(n_days_rough), time_years = day / 252, log_vol = logvol_path, vol = vol_path, H_true = H_true_rough, rough_logvol_scale = rough_logvol_scale)
rough_lags <- c(1, 2, 5, 10, 21, 42, 63)
rough_lags <- rough_lags[rough_lags < n_days_rough / 2]
rough_moment_tbl <- bind_rows(lapply(rough_lags, function(L) {
  diffs <- abs(logvol_path[(L + 1L):n_days_rough] - logvol_path[1:(n_days_rough - L)])
  tibble(lag_days = L, delta = L / 252, m1 = mean(diffs), m2 = mean(diffs^2))
}))
fit_m1 <- lm(log(m1) ~ log(delta), data = rough_moment_tbl)
fit_m2 <- lm(log(m2) ~ log(delta), data = rough_moment_tbl)
H_hat_m1 <- unname(coef(fit_m1)[2])
H_hat_m2 <- unname(coef(fit_m2)[2]) / 2
rough_empirical_summary <- tibble(H_true = H_true_rough, H_hat_from_q1_moment = H_hat_m1, H_hat_from_q2_moment = H_hat_m2, n_days = n_days_rough)

write_csv(rough_path_tbl, file.path(out_dir, "simulation_v2_simulated_logvol_path.csv"))
write_csv(rough_moment_tbl, file.path(out_dir, "simulation_v2_logvol_roughness_moments.csv"))
write_csv(rough_empirical_summary, file.path(out_dir, "simulation_v2_logvol_roughness_summary.csv"))
print(rough_empirical_summary)

# ============================================================
# 7. Experiment 2: iVi Volterra-Heston Monte Carlo simulation
# ============================================================

# Michael--Schucany--Haas sampler for the inverse Gaussian increments.
rinvgauss_ms <- function(mu, lambda) {
  mu <- pmax(mu, 1e-14)
  lambda <- pmax(lambda, 1e-14)
  y <- rnorm(length(mu))^2
  x <- mu + (mu^2 * y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
  u <- runif(length(mu))
  ifelse(u <= mu / (mu + x), x, mu^2 / x)
}

# iVi simulation of integrated variance increments and terminal stock values.
ivi_terminal_ST <- function(par, T, F = 100, n_steps = 64, n_paths = 20000) {
  alpha <- par$alpha
  dt <- T / n_steps
  b <- -par$gamma
  ccoef <- par$gamma * par$nu
  rho <- par$rho
  ell <- 0:(n_steps - 1L)
  k_int <- (((ell + 1) * dt)^alpha - (ell * dt)^alpha) / base::gamma(alpha + 1)
  k0 <- k_int[1L]
  Uinc <- matrix(0, nrow = n_paths, ncol = n_steps)
  Zinc <- matrix(0, nrow = n_paths, ncol = n_steps)
  alpha_min <- Inf
  alpha_clip_count <- 0L
  for (i in 0:(n_steps - 1L)) {
    ti <- i * dt
    tip1 <- (i + 1L) * dt
    alpha_i <- rep(par$V0 * dt + par$gamma * par$theta * (tip1^(alpha + 1) - ti^(alpha + 1)) / base::gamma(alpha + 2), n_paths)
    if (i > 0L) {
      for (j in 0:(i - 1L)) {
        lag <- i - j
        alpha_i <- alpha_i + k_int[lag + 1L] * (b * Uinc[, j + 1L] + ccoef * Zinc[, j + 1L])
      }
    }
    alpha_min <- min(alpha_min, min(alpha_i, na.rm = TRUE))
    bad_alpha <- !is.finite(alpha_i) | alpha_i <= 0
    alpha_clip_count <- alpha_clip_count + sum(bad_alpha)
    alpha_i <- pmax(alpha_i, 1e-14)
    mu_ig <- alpha_i / (1 - b * k0)
    lambda_ig <- (alpha_i / (ccoef * k0))^2
    Uinc[, i + 1L] <- rinvgauss_ms(mu_ig, lambda_ig)
    Zinc[, i + 1L] <- ((1 - b * k0) * Uinc[, i + 1L] - alpha_i) / (ccoef * k0)
  }
  if (alpha_clip_count > 0L) warning("iVi alpha_i clipping occurred. Inspect diagnostics before using results.")
  Nind <- matrix(rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
  dM_stock <- rho * Zinc + sqrt(1 - rho^2) * sqrt(pmax(Uinc, 0)) * Nind
  logS_T <- log(F) + rowSums(-0.5 * Uinc + dM_stock)
  ST <- exp(logS_T)
  diagnostics <- tibble(
    T = T,
    tau_days = round(365 * T),
    n_steps = n_steps,
    n_paths = n_paths,
    alpha_min_before_clip = alpha_min,
    alpha_clip_count = alpha_clip_count,
    alpha_clip_fraction = alpha_clip_count / (n_paths * n_steps),
    mean_ST = mean(ST),
    martingale_ratio_mean_ST_over_F = mean(ST) / F,
    martingale_error = mean(ST) / F - 1,
    sd_ST = sd(ST)
  )
  list(ST = ST, diagnostics = diagnostics)
}

run_ivi_surface_once <- function(par, Tval, tau_days, n_steps, n_paths, k_grid, F = F0, D = D0) {
  res <- ivi_terminal_ST(par = par, T = Tval, F = F, n_steps = n_steps, n_paths = n_paths)
  ST <- res$ST
  surface <- bind_rows(lapply(k_grid, function(kval) {
    K <- F * exp(kval)
    payoff <- pmax(ST - K, 0)
    mc_price <- D * mean(payoff)
    mc_price_se <- D * sd(payoff) / sqrt(length(payoff))
    mc_iv <- implied_vol_black_call(mc_price, F, K, Tval, D)
    vega_at_iv <- ifelse(is.finite(mc_iv), black_vega_forward(F, K, Tval, mc_iv, D), NA_real_)
    mc_iv_se <- ifelse(is.finite(vega_at_iv) && vega_at_iv > 1e-12, mc_price_se / vega_at_iv, NA_real_)
    tibble(tau = Tval, tau_days = tau_days, k = kval, K = K, n_steps = n_steps, n_paths = n_paths, mc_price = mc_price, mc_price_se = mc_price_se, mc_iv = mc_iv, mc_iv_se = mc_iv_se, H = par$H, alpha = par$alpha)
  }))
  list(surface = surface, diagnostics = res$diagnostics)
}

if (RUN_IVI_MODEL_EXPERIMENT) {
  cat("Running dedicated iVi Volterra-Heston Monte Carlo experiment...\n")
  ivi_base <- parameter_scenarios %>% filter(scenario_id == "iVi_stable_base") %>% slice(1)
  ivi_par <- list(gamma = ivi_base$gamma, theta = ivi_base$theta, nu = ivi_base$nu, rho = ivi_base$rho, V0 = ivi_base$V0, alpha = ivi_base$alpha)
  ivi_par$H <- ivi_par$alpha - 0.5
  ivi_parameter_tbl <- tibble(
    scenario_id = "iVi_stable_base",
    gamma = ivi_par$gamma,
    theta = ivi_par$theta,
    nu = ivi_par$nu,
    rho = ivi_par$rho,
    V0 = ivi_par$V0,
    alpha = ivi_par$alpha,
    H = ivi_par$H,
    n_paths_main = IVI_N_PATHS,
    n_paths_reference = IVI_REF_N_PATHS,
    reference_n_steps = IVI_REF_N_STEPS,
    tau_days_grid = paste(IVI_TAU_DAYS, collapse = ", "),
    k_grid = paste(IVI_K_GRID, collapse = ", "),
    parameter_status = "model_based_simulation_benchmark_not_market_calibrated",
    fast_test = FAST_TEST
  )
  write_csv(ivi_parameter_tbl, file.path(out_dir, "simulation_v2_ivi_parameters.csv"))
  print(ivi_parameter_tbl)
  ivi_runs <- list(); ivi_diag_runs <- list(); idx <- 1L
  for (td in IVI_TAU_DAYS) {
    Tval <- td / 365
    for (ns in IVI_N_STEPS_GRID) {
      cat("  iVi main run: tau_days =", td, "| n_steps =", ns, "\n")
      one <- run_ivi_surface_once(ivi_par, Tval, td, ns, IVI_N_PATHS, IVI_K_GRID, F0, D0)
      ivi_runs[[idx]] <- one$surface
      ivi_diag_runs[[idx]] <- one$diagnostics
      idx <- idx + 1L
    }
  }
  ivi_mc_surface <- bind_rows(ivi_runs)
  ivi_mc_diagnostics <- bind_rows(ivi_diag_runs)
  ivi_ref_runs <- list(); ivi_ref_diag_runs <- list(); ridx <- 1L
  for (td in IVI_TAU_DAYS) {
    Tval <- td / 365
    cat("  iVi reference run: tau_days =", td, "| n_steps =", IVI_REF_N_STEPS, "\n")
    one_ref <- run_ivi_surface_once(ivi_par, Tval, td, IVI_REF_N_STEPS, IVI_REF_N_PATHS, IVI_K_GRID, F0, D0)
    ivi_ref_runs[[ridx]] <- one_ref$surface
    ivi_ref_diag_runs[[ridx]] <- one_ref$diagnostics
    ridx <- ridx + 1L
  }
  ivi_reference <- bind_rows(ivi_ref_runs)
  ivi_reference_diagnostics <- bind_rows(ivi_ref_diag_runs)
  ivi_surface_comparison <- ivi_mc_surface %>%
    left_join(
      ivi_reference %>% select(tau_days, k, ref_n_steps = n_steps, ref_n_paths = n_paths, ref_price = mc_price, ref_price_se = mc_price_se, ref_iv = mc_iv, ref_iv_se = mc_iv_se),
      by = c("tau_days", "k")
    ) %>%
    mutate(iv_error_to_ref = mc_iv - ref_iv, abs_iv_error_to_ref = abs(iv_error_to_ref), combined_iv_se = sqrt(mc_iv_se^2 + ref_iv_se^2))
  ivi_convergence_summary <- ivi_surface_comparison %>%
    group_by(tau_days, tau, n_steps, n_paths) %>%
    summarise(
      n_quotes = sum(is.finite(iv_error_to_ref)),
      iv_rmse_to_reference = sqrt(mean(iv_error_to_ref^2, na.rm = TRUE)),
      iv_mae_to_reference = mean(abs(iv_error_to_ref), na.rm = TRUE),
      mean_mc_iv_se = mean(mc_iv_se, na.rm = TRUE),
      mean_ref_iv_se = mean(ref_iv_se, na.rm = TRUE),
      mean_combined_iv_se = mean(combined_iv_se, na.rm = TRUE),
      max_abs_iv_error_to_reference = max(abs_iv_error_to_ref, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      rmse_over_combined_mcse = iv_rmse_to_reference / mean_combined_iv_se,
      convergence_comment = case_when(
        rmse_over_combined_mcse <= 1.0 ~ "discretisation error below combined MC noise",
        rmse_over_combined_mcse <= 1.5 ~ "discretisation error near combined MC noise",
        TRUE ~ "discretisation error visible or MC reference still noisy"
      )
    )
  ivi_ref_atm_skew <- ivi_reference %>%
    group_by(tau_days, tau) %>%
    group_modify(~ estimate_atm_skew(.x, "mc_iv", window = 0.055)) %>%
    ungroup() %>%
    rename(
      ivi_ref_atm_skew_lq = psi,
      ivi_ref_atm_iv_lq = intercept,
      ivi_ref_n_local = n_local,
      ivi_ref_fit_type = fit_type
    )
  
  ivi_ref_support <- ivi_reference %>%
    group_by(tau_days, tau) %>%
    group_modify(~ near_atm_support_audit(.x, window = 0.055)) %>%
    ungroup()
  
  ivi_ref_cd_h0025 <- ivi_reference %>%
    group_by(tau_days, tau) %>%
    group_modify(~ central_difference_skew(.x, "mc_iv", h = 0.025)) %>%
    ungroup() %>%
    rename(ivi_ref_cd_skew_h0025 = cd_skew, ivi_ref_cd_available_h0025 = cd_available)
  
  ivi_ref_cd_h005 <- ivi_reference %>%
    group_by(tau_days, tau) %>%
    group_modify(~ central_difference_skew(.x, "mc_iv", h = 0.05)) %>%
    ungroup() %>%
    rename(ivi_ref_cd_skew_h005 = cd_skew, ivi_ref_cd_available_h005 = cd_available)
  
  ivi_ref_skew_audit <- ivi_ref_atm_skew %>%
    left_join(ivi_ref_support, by = c("tau_days", "tau")) %>%
    left_join(ivi_ref_cd_h0025, by = c("tau_days", "tau")) %>%
    left_join(ivi_ref_cd_h005, by = c("tau_days", "tau")) %>%
    mutate(
      headline_skew_estimator_ok = two_sided_support & ivi_ref_fit_type == "local_quadratic" & ivi_ref_n_local >= 5L,
      lq_minus_cd_h0025 = ivi_ref_atm_skew_lq - ivi_ref_cd_skew_h0025,
      lq_minus_cd_h005 = ivi_ref_atm_skew_lq - ivi_ref_cd_skew_h005
    )
  
  write_csv(ivi_mc_surface, file.path(out_dir, "simulation_v2_ivi_mc_surface.csv"))
  write_csv(ivi_reference, file.path(out_dir, "simulation_v2_ivi_mc_reference_surface.csv"))
  write_csv(ivi_surface_comparison, file.path(out_dir, "simulation_v2_ivi_mc_surface_comparison.csv"))
  write_csv(ivi_convergence_summary, file.path(out_dir, "simulation_v2_ivi_mc_convergence_summary.csv"))
  write_csv(ivi_convergence_summary, file.path(out_dir, "simulation_v2_ivi_mc_convergence_diagnostic.csv"))
  write_csv(ivi_mc_diagnostics, file.path(out_dir, "simulation_v2_ivi_mc_diagnostics.csv"))
  write_csv(ivi_reference_diagnostics, file.path(out_dir, "simulation_v2_ivi_mc_reference_diagnostics.csv"))
  write_csv(ivi_ref_atm_skew, file.path(out_dir, "simulation_v2_ivi_reference_atm_skew_lq.csv"))
  write_csv(ivi_ref_skew_audit, file.path(out_dir, "simulation_v2_ivi_reference_atm_skew_extraction_audit.csv"))
  print(ivi_convergence_summary)
  print(ivi_reference_diagnostics)
  print(ivi_ref_skew_audit)
}

# ============================================================
# 8. Plots
# ============================================================

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5, colour = "grey20"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
    plot.margin = margin(8, 10, 8, 10)
  )

maturity_palette <- c("7d"="#1B9E77","14d"="#D95F02","21d"="#7570B3","30d"="#E7298A","45d"="#66A61E","60d"="#E6AB02","90d"="#A6761D","120d"="#1F78B4","180d"="#B2DF8A","270d"="#FB9A99","365d"="#6A3D9A")

smile_days <- intersect(c(7, 14, 30, 60, 90, 180, 365), unique(surface_grid$tau_days))
p_smiles <- surface_grid %>%
  filter(tau_days %in% smile_days) %>%
  mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(smile_days, "d"))) %>%
  arrange(tau_days, k) %>%
  ggplot(aes(x = k, y = target_iv, colour = tau_label)) +
  geom_line(linewidth = 0.75, show.legend = FALSE) +
  geom_point(size = 1.05, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey35") +
  facet_wrap(~ tau_label, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = maturity_palette) +
  labs(title = "Synthetic implied volatility smiles", x = "Forward log-moneyness, k = log(K/F)", y = "B&S implied volatility") + 
  plot_theme
plot_save(p_smiles, "plot_v2_synthetic_smiles", width = 11, height = 7)

loglog_line <- tibble(tau = seq(min(atm_skew_tbl$tau), max(atm_skew_tbl$tau), length.out = 200)) %>% mutate(abs_psi_fit = exp(atm_intercept) * tau^atm_slope)
p_loglog <- atm_skew_tbl %>%
  filter(abs_psi > 0, is.finite(abs_psi)) %>%
  ggplot(aes(x = tau, y = abs_psi)) +
  geom_point(size = 2.2, colour = "#0072B2") +
  geom_line(data = loglog_line, aes(x = tau, y = abs_psi_fit), inherit.aes = FALSE, linewidth = 0.75, colour = "#D55E00") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Log-log scaling of the ATM skew", x = expression(tau), y = expression(abs(psi(tau)))) +
  plot_theme
plot_save(p_loglog, "plot_v2_loglog_atm_skew", width = 8.5, height = 5.2)

if (RUN_IVI_MODEL_EXPERIMENT && exists("ivi_mc_surface") && exists("ivi_reference")) {
  ivi_main_plot_tbl <- ivi_mc_surface %>% mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(IVI_TAU_DAYS, "d")), step_label = factor(paste0("n=", n_steps), levels = paste0("n=", IVI_N_STEPS_GRID))) %>% arrange(tau_days, n_steps, k)
  ivi_ref_plot_tbl <- ivi_reference %>% mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(IVI_TAU_DAYS, "d"))) %>% arrange(tau_days, k)
  p_ivi_smiles <- ggplot() +
    geom_line(data = ivi_main_plot_tbl, aes(x = k, y = mc_iv, colour = step_label, group = step_label), linewidth = 0.45, alpha = 0.85) +
    geom_point(data = ivi_main_plot_tbl, aes(x = k, y = mc_iv, colour = step_label), size = 0.85, alpha = 0.85) +
    geom_line(data = ivi_ref_plot_tbl, aes(x = k, y = mc_iv), linewidth = 0.95, colour = "#000000") +
    geom_point(data = ivi_ref_plot_tbl, aes(x = k, y = mc_iv), size = 1.25, colour = "#000000") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey35") +
    facet_wrap(~ tau_label, scales = "free_y") +
    labs(title = "iVi Volterra-Heston implied volatility slices", x = "Forward log-moneyness, k = log(K/F)", y = "Monte Carlo B&S implied volatility", colour = "iVi time steps") +
    plot_theme
  plot_save(p_ivi_smiles, "plot_v2_ivi_mc_smiles", width = 11, height = 6.2)
  p_ivi_convergence <- ivi_convergence_summary %>% mutate(tau_label = factor(paste0(tau_days, "d"), levels = paste0(IVI_TAU_DAYS, "d"))) %>%
    ggplot(aes(x = n_steps, y = iv_rmse_to_reference, colour = tau_label, linetype = tau_label)) +
    geom_line(linewidth = 0.75) + geom_point(size = 2.0) + scale_x_log10(breaks = IVI_N_STEPS_GRID) +
    labs(title = "iVi time-step convergence", x = "Number of time steps", y = "IV RMSE to reference", colour = "Maturity", linetype = "Maturity") +
    plot_theme
  plot_save(p_ivi_convergence, "plot_v2_ivi_mc_convergence", width = 8.5, height = 5.2)
}

# Figure manifest separates report figures from diagnostic-only plots.
figure_manifest <- tibble(
  file = paste0(c("plot_v2_synthetic_smiles", "plot_v2_loglog_atm_skew", "plot_v2_ivi_mc_smiles", "plot_v2_ivi_mc_convergence"), ".", PLOT_FORMAT),
  use_in_thesis = c(TRUE, TRUE, TRUE, FALSE),
  description = c(
    "Controlled SPX-scale synthetic IV smiles with negative equity skew.",
    "Log-log slope check for |psi(tau)| against H - 1/2.",
    "Model-based iVi Volterra-Heston Monte Carlo implied volatility slices.",
    "iVi time-step convergence diagnostic against the reference Monte Carlo run."
  )
)
write_csv(figure_manifest, file.path(out_dir, "simulation_v2_figure_manifest.csv"))

cat("\n=== Phase S completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("FAST_TEST mode:", FAST_TEST, "\n")
cat("Synthetic reduced-form surface quotes:", nrow(surface_grid), "\n")
cat("Calibration quotes:", nrow(calib_panel), "\n")
cat("Held-out test quotes:", nrow(test_panel), "\n")
cat("ATM skew slope estimate:", round(atm_slope, 4), "\n")
cat("Reference H - 1/2:", round(rough_design$H - 0.5, 4), "\n")
cat("Recovered H from ATM skew slope:", round(atm_summary$H_recovered_from_slope[1], 4), "\n")
cat("Recovered H from reduced-form full-surface fit:", round(rough_hat$H, 4), "\n")
cat("Simulated log-vol H estimate from q=1 moment:", round(H_hat_m1, 4), "\n")
cat("Simulated log-vol H estimate from q=2 moment:", round(H_hat_m2, 4), "\n")
cat("Held-out rough IV RMSE vs target:", signif(error_summary$rough_iv_rmse_target[1], 4), "\n")
cat("Held-out smooth H=1/2 IV RMSE vs target:", signif(error_summary$smooth_iv_rmse_target[1], 4), "\n")
if (RUN_IVI_MODEL_EXPERIMENT && exists("ivi_convergence_summary")) {
  cat("iVi model-based experiment: completed.\n")
  cat("iVi maturities:", paste(IVI_TAU_DAYS, collapse = ", "), "days\n")
  cat("iVi main paths:", IVI_N_PATHS, "\n")
  cat("iVi reference paths:", IVI_REF_N_PATHS, "\n")
}
cat("Files written: CSV tables and report-ready plots.\n")
cat("Next step: run Phase_s_rheston_cf_pricing.R.\n")
