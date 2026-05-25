#!/usr/bin/env Rscript

# ============================================================
# Empirical: 2010 SPX market panel construction
# ------------------------------------------------------------
#
# Purpose:
#   Build a strict one-date SPX empirical calibration panel for
#   the main rough Heston model and benchmark comparisons against
#   classical Heston and a pure-jump Variance Gamma model.
#
# Inputs:
#   spx_option_rawiv_full_2010_01_07.csv
#   spx_full_option_prices_2010_01_07.csv       
#   spx_stock_prices_2010_01_07.csv
#   usd_interest_rates_2010_01_07.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260107)

# ============================================================
# 0. Run controls and paths
# ============================================================

FAST_TEST <- identical(Sys.getenv("PHASEM_FAST_TEST"), "1")
CLEAN_OUTPUT_DIR <- TRUE
PLOT_FORMAT <- "pdf"
SHOW_PLOTS <- FALSE
SAVE_PLOTS <- TRUE

target_date <- as.Date("2010-01-07")
date_tag <- "2010_01_07"

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

candidate_dirs <- unique(c(
  file.path(base_dir, "data", "Benchmark"),
  file.path(base_dir, "data"),
  base_dir,
  getwd(),
  "/mnt/data"
))

find_input_file <- function(fname) {
  hits <- file.path(candidate_dirs, fname)
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0L) {
    stop("Could not find input file: ", fname,
         "\nSearched:\n", paste(candidate_dirs, collapse = "\n"))
  }
  hits[1]
}

out_dir <- file.path(base_dir, "data", "EmpiricalStudy", "phaseM_2010_market_panel_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_OUTPUT_DIR) {
  old_files <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
  if (length(old_files) > 0L) unlink(old_files, recursive = TRUE, force = TRUE)
}

cat("Empirical 2010 market-panel output directory:\n", out_dir, "\n\n", sep = "")
cat("FAST_TEST mode:", FAST_TEST, "\n\n")

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

# ============================================================
# 1. Black-Scholes functions on forwards
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

black_put_forward <- function(F, K, tau, sigma, D = 1) {
  if (!is.finite(sigma) || sigma <= 0 || tau <= 0) {
    return(D * max(K - F, 0))
  }
  vol_sqrt <- sigma * sqrt(tau)
  d1 <- (log(F / K) + 0.5 * sigma^2 * tau) / vol_sqrt
  d2 <- d1 - vol_sqrt
  D * (K * pnorm(-d2) - F * pnorm(-d1))
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
# 2. Read and normalize data
# ============================================================

rawiv_path <- find_input_file(paste0("spx_option_rawiv_full_", date_tag, ".csv"))
prices_path <- find_input_file(paste0("spx_full_option_prices_", date_tag, ".csv"))
stock_path <- find_input_file(paste0("spx_stock_prices_", date_tag, ".csv"))
rates_path <- find_input_file(paste0("usd_interest_rates_", date_tag, ".csv"))

rawiv <- read_csv(rawiv_path, show_col_types = FALSE)
prices_raw <- read_csv(prices_path, show_col_types = FALSE)
stock <- read_csv(stock_path, show_col_types = FALSE)
rates <- read_csv(rates_path, show_col_types = FALSE)

S0 <- as.numeric(stock$close[1])
if (!is.finite(S0) || S0 <= 0) {
  S0 <- as.numeric(rawiv$`Adjusted close`[1])
}

rates_clean <- rates %>%
  transmute(
    period_days = as.numeric(period),
    rate_decimal = as.numeric(rate) / 100
  ) %>%
  filter(is.finite(period_days), is.finite(rate_decimal)) %>%
  arrange(period_days) %>%
  distinct(period_days, .keep_all = TRUE)

rate_for_days <- function(days) {
  approx(rates_clean$period_days, rates_clean$rate_decimal, xout = days, rule = 2)$y
}

opts0 <- rawiv %>%
  rename(
    option_symbol = `option symbol`,
    cp = `Call/Put`,
    K = strike,
    adj_close = `Adjusted close`,
    raw_iv = iv,
    open_interest = `open interest`
  ) %>%
  mutate(
    date = as.Date(date),
    expiration = as.Date(expiration),
    
    # Keep the vendor-listed expiration, but use the preceding Friday
    # as the pricing maturity when standard SPX monthly expiry is stored
    # as a Saturday. POSIX weekday convention: 0 = Sunday, ..., 6 = Saturday.
    expiration_wday = as.POSIXlt(expiration)$wday,
    pricing_expiration = ifelse(expiration_wday == 6L, expiration - 1, expiration),
    pricing_expiration = as.Date(pricing_expiration, origin = "1970-01-01"),
    
    tau_days = as.integer(pricing_expiration - date),
    tau = tau_days / 365,
    cp = toupper(cp),
    K = as.numeric(K),
    bid = as.numeric(bid),
    ask = as.numeric(ask),
    price_vendor = as.numeric(price),
    mid = 0.5 * (bid + ask),
    spread = ask - bid,
    rel_spread = spread / pmax(mid, 1e-8),
    raw_iv = as.numeric(raw_iv),
    volume = as.numeric(volume),
    open_interest = as.numeric(open_interest)
  ) %>%
  filter(
    date == target_date,
    cp %in% c("C", "P"),
    tau_days > 0,
    is.finite(K), K > 0,
    is.finite(bid), is.finite(ask),
    ask >= bid,
    ask > 0,
    is.finite(mid), mid > 0,
    is.finite(raw_iv), raw_iv > 0
  )

raw_audit <- tibble(
  valuation_date = target_date,
  S0 = S0,
  raw_rows = nrow(rawiv),
  normalized_valid_rows = nrow(opts0),
  n_expiries = n_distinct(opts0$expiration),
  min_tau_days = min(opts0$tau_days),
  max_tau_days = max(opts0$tau_days),
  n_calls = sum(opts0$cp == "C"),
  n_puts = sum(opts0$cp == "P")
)

write_csv(raw_audit, file.path(out_dir, "phaseM2010_raw_audit.csv"))
print(raw_audit)

# ============================================================
# 3. Put-call parity forward and discount extraction
# ============================================================

opts_one <- opts0 %>%
  group_by(expiration, tau_days, tau, K, cp) %>%
  summarise(
    mid = mean(mid, na.rm = TRUE),
    bid = mean(bid, na.rm = TRUE),
    ask = mean(ask, na.rm = TRUE),
    spread = mean(spread, na.rm = TRUE),
    raw_iv = mean(raw_iv, na.rm = TRUE),
    volume = sum(volume, na.rm = TRUE),
    open_interest = sum(open_interest, na.rm = TRUE),
    .groups = "drop"
  )

pairs <- opts_one %>%
  select(expiration, tau_days, tau, K, cp, mid, bid, ask, spread, raw_iv, volume, open_interest) %>%
  pivot_wider(
    names_from = cp,
    values_from = c(mid, bid, ask, spread, raw_iv, volume, open_interest),
    names_sep = "_"
  ) %>%
  filter(is.finite(mid_C), is.finite(mid_P))

estimate_parity_one <- function(df, S0) {
  d0 <- df %>%
    mutate(
      y = mid_C - mid_P,
      spread_sum = pmax(spread_C + spread_P, 0.01),
      log_moneyness_spot = log(K / S0)
    ) %>%
    filter(is.finite(y), is.finite(K), is.finite(spread_sum))
  
  # Use central strikes for parity; widen if necessary.
  d <- d0 %>% filter(abs(log_moneyness_spot) <= 0.30)
  if (nrow(d) < 8L) d <- d0 %>% filter(abs(log_moneyness_spot) <= 0.50)
  if (nrow(d) < 5L) {
    return(tibble(
      n_pairs = nrow(d0),
      n_pairs_used = nrow(d),
      F_parity = NA_real_,
      D_parity = NA_real_,
      parity_rmse = NA_real_,
      parity_ok = FALSE
    ))
  }
  
  fit <- tryCatch(
    lm(y ~ K, data = d, weights = 1 / pmax(spread_sum^2, 1e-4)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      n_pairs = nrow(d0),
      n_pairs_used = nrow(d),
      F_parity = NA_real_,
      D_parity = NA_real_,
      parity_rmse = NA_real_,
      parity_ok = FALSE
    ))
  }
  
  cf <- coef(fit)
  intercept <- unname(cf["(Intercept)"])
  slope <- unname(cf["K"])
  D_hat <- -slope
  F_hat <- intercept / D_hat
  resid <- residuals(fit)
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  
  parity_ok <- is.finite(F_hat) && is.finite(D_hat) &&
    F_hat > 0.5 * S0 && F_hat < 1.5 * S0 &&
    D_hat > 0.80 && D_hat < 1.05 &&
    is.finite(rmse)
  
  tibble(
    n_pairs = nrow(d0),
    n_pairs_used = nrow(d),
    F_parity = F_hat,
    D_parity = D_hat,
    parity_rmse = rmse,
    parity_ok = parity_ok
  )
}

parity_tbl <- pairs %>%
  group_by(expiration, tau_days, tau) %>%
  group_modify(~ estimate_parity_one(.x, S0 = S0)) %>%
  ungroup() %>%
  mutate(
    r_curve = rate_for_days(tau_days),
    D_curve = exp(-r_curve * tau),
    F_curve_no_div = S0 / D_curve,
    D_abs_diff_curve = abs(D_parity - D_curve),
    D_rel_diff_curve = D_abs_diff_curve / pmax(D_curve, 1e-12)
  ) %>%
  arrange(tau_days)

write_csv(parity_tbl, file.path(out_dir, "phaseM2010_put_call_parity_by_expiry.csv"))
print(parity_tbl)

if (any(!parity_tbl$parity_ok, na.rm = FALSE)) {
  warning("At least one expiry failed parity extraction. Inspect phaseM2010_put_call_parity_by_expiry.csv.")
}

# ============================================================
# 4. Build OTM call-equivalent market panel
# ============================================================

opts_panel0 <- opts0 %>%
  left_join(
    parity_tbl %>%
      select(expiration, tau_days, tau, F = F_parity, D = D_parity, parity_ok, parity_rmse),
    by = c("expiration", "tau_days", "tau")
  ) %>%
  filter(parity_ok, is.finite(F), is.finite(D), F > 0, D > 0) %>%
  mutate(
    k = log(K / F),
    is_otm = (cp == "C" & K >= F) | (cp == "P" & K < F),
    otm_score = ifelse(is_otm, 0, 1) + 0.01 * pmin(rel_spread, 100)
  ) %>%
  group_by(expiration, tau_days, tau, K) %>%
  arrange(otm_score, rel_spread, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(
    selected_is_otm = is_otm,
    call_equiv_price = ifelse(cp == "C", mid, mid + D * (F - K)),
    call_equiv_bid = ifelse(cp == "C", bid, bid + D * (F - K)),
    call_equiv_ask = ifelse(cp == "C", ask, ask + D * (F - K)),
    intrinsic_call = D * pmax(F - K, 0),
    upper_call = D * F,
    call_equiv_valid = is.finite(call_equiv_price) &
      call_equiv_price >= intrinsic_call - 1e-8 &
      call_equiv_price <= upper_call + 1e-8
  )

# Recompute the calibration implied volatility from the call-equivalent mid price.
# The vendor raw IV is retained as a diagnostic field only; it is not used
# as a fallback calibration target.
iv_calc <- mapply(
  function(price, F, K, tau, D) implied_vol_black_call(price, F, K, tau, D),
  opts_panel0$call_equiv_price,
  opts_panel0$F,
  opts_panel0$K,
  opts_panel0$tau,
  opts_panel0$D
)

opts_panel <- opts_panel0 %>%
  mutate(
    iv_calc = as.numeric(iv_calc),
    iv_mkt = iv_calc,
    vega_mkt = mapply(
      function(F_, K_, tau_, iv_, D_) black_vega_forward(F_, K_, tau_, iv_, D_),
      F, K, tau, iv_mkt, D
    ),
    iv_spread_est = spread / pmax(vega_mkt, 1e-8),
    iv_spread_est = pmin(pmax(iv_spread_est, 0.0005), 2.0)
  ) %>%
  filter(
    call_equiv_valid,
    is.finite(iv_mkt), iv_mkt > 0.01, iv_mkt < 3.0,
    is.finite(k), abs(k) <= 0.80
  )

# Diagnostic only: compare recomputed Black--Scholes implied volatility
# against the vendor raw IV on the cleaned calibration convention.
iv_recompute_diagnostic <- opts_panel %>%
  mutate(iv_abs_diff = abs(iv_calc - raw_iv)) %>%
  group_by(expiration, tau_days, tau) %>%
  summarise(
    n_quotes = n(),
    median_abs_diff = median(iv_abs_diff, na.rm = TRUE),
    p95_abs_diff = as.numeric(quantile(iv_abs_diff, 0.95, na.rm = TRUE)),
    max_abs_diff = max(iv_abs_diff, na.rm = TRUE),
    share_abs_diff_below_1volbp = mean(iv_abs_diff <= 0.0001, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tau_days)
write_csv(iv_recompute_diagnostic, file.path(out_dir, "phaseM2010_rawiv_recomputed_iv_diagnostic.csv"))

# Weight construction. Trim extremely large weights.
w_raw <- 1 / pmax(opts_panel$iv_spread_est^2, 1e-6)
w_cap <- as.numeric(quantile(w_raw[is.finite(w_raw)], 0.95, na.rm = TRUE))
opts_panel <- opts_panel %>%
  mutate(
    weight_iv = pmin(w_raw, w_cap),
    tau_label = paste0(tau_days, "d")
  )

# ============================================================
# 5. ATM-skew extraction, support audit, and H selection
# ============================================================

estimate_atm_skew <- function(df, iv_col = "iv_mkt", window = 0.055) {
  d <- df %>%
    filter(abs(k) <= window, is.finite(.data[[iv_col]])) %>%
    arrange(k)
  
  n_left <- sum(d$k < -1e-10)
  n_atm <- sum(abs(d$k) <= 0.002)
  n_right <- sum(d$k > 1e-10)
  two_sided_support <- n_left >= 1L && n_right >= 1L
  k_support <- if (nrow(d) > 0L) paste(sprintf("%.4f", d$k), collapse = ", ") else ""
  
  if (nrow(d) < 3L || !two_sided_support) {
    return(tibble(
      atm_iv = NA_real_,
      atm_skew = NA_real_,
      n_local = nrow(d),
      n_left = n_left,
      n_atm = n_atm,
      n_right = n_right,
      two_sided_support = two_sided_support,
      fit_type = "insufficient",
      k_support = k_support,
      skew_quality = "reject_insufficient_two_sided_support"
    ))
  }
  
  d <- d %>% mutate(local_weight = pmax(0, 1 - abs(k) / window))
  
  # Match the validated simulation hierarchy: local quadratic is the
  # headline estimator when the near-ATM support is sufficiently rich and
  # two-sided. Local linear estimates are retained only as lower-quality
  # diagnostics and are not used for the headline H estimate.
  fit_type <- if (nrow(d) >= 5L && n_left >= 2L && n_right >= 2L) {
    "local_quadratic"
  } else {
    "local_linear"
  }
  
  form <- if (fit_type == "local_quadratic") {
    as.formula(paste0(iv_col, " ~ k + I(k^2)"))
  } else {
    as.formula(paste0(iv_col, " ~ k"))
  }
  
  fit <- tryCatch(lm(form, data = d, weights = local_weight), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      atm_iv = NA_real_,
      atm_skew = NA_real_,
      n_local = nrow(d),
      n_left = n_left,
      n_atm = n_atm,
      n_right = n_right,
      two_sided_support = two_sided_support,
      fit_type = "failed",
      k_support = k_support,
      skew_quality = "reject_fit_failed"
    ))
  }
  
  cf <- coef(fit)
  quality <- if (fit_type == "local_quadratic") {
    "headline_local_quadratic"
  } else {
    "lower_quality_local_linear"
  }
  
  tibble(
    atm_iv = unname(cf["(Intercept)"]),
    atm_skew = unname(cf["k"]),
    n_local = nrow(d),
    n_left = n_left,
    n_atm = n_atm,
    n_right = n_right,
    two_sided_support = two_sided_support,
    fit_type = fit_type,
    k_support = k_support,
    skew_quality = quality
  )
}

central_difference_one <- function(df, iv_col = "iv_mkt", h = 0.025, max_gap = 0.020) {
  d <- df %>% filter(is.finite(k), is.finite(.data[[iv_col]]))
  
  # In market data the strike grid is not exactly symmetric in forward
  # log-moneyness. This robustness check therefore uses the closest
  # available negative and positive points to -h and +h, provided that both
  # are sufficiently close to the requested locations.
  kp <- d %>% filter(k > 0) %>% slice_min(abs(k - h), n = 1, with_ties = FALSE)
  km <- d %>% filter(k < 0) %>% slice_min(abs(k + h), n = 1, with_ties = FALSE)
  
  if (nrow(kp) == 0L || nrow(km) == 0L || abs(kp$k[1] - h) > max_gap || abs(km$k[1] + h) > max_gap) {
    return(tibble(
      h = h,
      cd_available = FALSE,
      cd_skew = NA_real_,
      k_minus = ifelse(nrow(km) > 0L, km$k[1], NA_real_),
      k_plus = ifelse(nrow(kp) > 0L, kp$k[1], NA_real_),
      h_effective = NA_real_
    ))
  }
  
  tibble(
    h = h,
    cd_available = TRUE,
    cd_skew = (kp[[iv_col]][1] - km[[iv_col]][1]) / (kp$k[1] - km$k[1]),
    k_minus = km$k[1],
    k_plus = kp$k[1],
    h_effective = 0.5 * (kp$k[1] - km$k[1])
  )
}

sign_same <- function(x, y) {
  is.finite(x) & is.finite(y) & sign(x) == sign(y) & sign(x) != 0
}

fit_loglog_skew <- function(df, label, H_floor = 0.03, H_cap = 0.49) {
  d <- df %>%
    filter(is.finite(atm_skew), abs(atm_skew) > 0, is.finite(tau), tau > 0) %>%
    mutate(abs_atm_skew = abs(atm_skew))
  
  if (nrow(d) < 3L) {
    return(list(
      summary = tibble(
        estimator = label,
        n_maturities = nrow(d),
        maturity_days_used = paste(d$tau_days, collapse = ", "),
        intercept = NA_real_, beta_hat = NA_real_, H_hat = NA_real_,
        H_for_calibration = NA_real_, alpha_for_calibration = NA_real_,
        H_was_clipped = NA, r_squared = NA_real_
      ),
      points = d %>% mutate(
        estimator = label,
        log_abs_atm_skew = log(abs_atm_skew),
        fitted_log_abs_atm_skew = NA_real_,
        loglog_residual = NA_real_
      )
    ))
  }
  
  fit <- lm(log(abs_atm_skew) ~ log(tau), data = d)
  beta_hat <- unname(coef(fit)[2])
  H_hat <- beta_hat + 0.5
  H_for_cal <- min(max(H_hat, H_floor), H_cap)
  
  points <- d %>%
    mutate(
      estimator = label,
      log_abs_atm_skew = log(abs_atm_skew),
      fitted_log_abs_atm_skew = fitted(fit),
      loglog_residual = residuals(fit)
    )
  
  summary <- tibble(
    estimator = label,
    n_maturities = nrow(d),
    maturity_days_used = paste(d$tau_days, collapse = ", "),
    intercept = unname(coef(fit)[1]),
    beta_hat = beta_hat,
    H_hat = H_hat,
    H_for_calibration = H_for_cal,
    alpha_for_calibration = H_for_cal + 0.5,
    H_was_clipped = abs(H_for_cal - H_hat) > 1e-12,
    r_squared = summary(fit)$r.squared
  )
  
  list(summary = summary, points = points)
}

atm_skew_tbl <- opts_panel %>%
  group_by(expiration, tau_days, tau) %>%
  group_modify(~ estimate_atm_skew(.x, "iv_mkt", window = 0.055)) %>%
  ungroup() %>%
  mutate(
    good_atm_skew = is.finite(atm_skew) & is.finite(atm_iv) & two_sided_support & n_local >= 3,
    headline_atm_skew = good_atm_skew & fit_type == "local_quadratic" & n_local >= 5L,
    very_short_maturity = tau_days <= 2L,
    very_short_support_decision = case_when(
      !very_short_maturity ~ "not_very_short",
      good_atm_skew ~ "passes_support_audit",
      TRUE ~ "fails_support_audit"
    )
  ) %>%
  arrange(tau_days)

near_atm_support_audit <- opts_panel %>%
  group_by(expiration, tau_days, tau) %>%
  summarise(
    n_panel_quotes = n(),
    n_near_atm = sum(abs(k) <= 0.055, na.rm = TRUE),
    n_left = sum(k < -1e-10 & abs(k) <= 0.055, na.rm = TRUE),
    n_atm_band = sum(abs(k) <= 0.002, na.rm = TRUE),
    n_right = sum(k > 1e-10 & abs(k) <= 0.055, na.rm = TRUE),
    k_min_near = ifelse(n_near_atm > 0, min(k[abs(k) <= 0.055], na.rm = TRUE), NA_real_),
    k_max_near = ifelse(n_near_atm > 0, max(k[abs(k) <= 0.055], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  left_join(atm_skew_tbl, by = c("expiration", "tau_days", "tau")) %>%
  arrange(tau_days)

cd_h0025 <- opts_panel %>%
  group_by(expiration, tau_days, tau) %>%
  group_modify(~ central_difference_one(.x, "iv_mkt", h = 0.025, max_gap = 0.020)) %>%
  ungroup() %>%
  rename(cd_available_h0025 = cd_available, cd_skew_h0025 = cd_skew, k_minus_h0025 = k_minus, k_plus_h0025 = k_plus, h_effective_h0025 = h_effective) %>%
  select(-h)

cd_h005 <- opts_panel %>%
  group_by(expiration, tau_days, tau) %>%
  group_modify(~ central_difference_one(.x, "iv_mkt", h = 0.05, max_gap = 0.020)) %>%
  ungroup() %>%
  rename(cd_available_h005 = cd_available, cd_skew_h005 = cd_skew, k_minus_h005 = k_minus, k_plus_h005 = k_plus, h_effective_h005 = h_effective) %>%
  select(-h)

atm_skew_cd_tbl <- atm_skew_tbl %>%
  select(expiration, tau_days, tau, atm_iv, atm_skew, n_local, fit_type, good_atm_skew, headline_atm_skew, very_short_maturity, very_short_support_decision) %>%
  left_join(cd_h0025, by = c("expiration", "tau_days", "tau")) %>%
  left_join(cd_h005, by = c("expiration", "tau_days", "tau")) %>%
  mutate(
    lq_minus_cd_h0025 = atm_skew - cd_skew_h0025,
    lq_minus_cd_h005 = atm_skew - cd_skew_h005,
    cd_sign_same_h0025 = sign_same(atm_skew, cd_skew_h0025),
    cd_sign_same_h005 = sign_same(atm_skew, cd_skew_h005),
    cd_robust_for_headline = headline_atm_skew & cd_available_h0025 & cd_available_h005 & cd_sign_same_h0025 & cd_sign_same_h005,
    cd_robust_comment = case_when(
      !headline_atm_skew ~ "not_headline_local_quadratic",
      !cd_available_h0025 | !cd_available_h005 ~ "central_difference_not_available",
      !cd_sign_same_h0025 | !cd_sign_same_h005 ~ "central_difference_sign_instability",
      TRUE ~ "central_difference_signs_consistent"
    )
  ) %>%
  arrange(tau_days)

# H selection for the roughness input:
#   * estimate H from the short-end ATM-skew term structure;
#   * audit the one-day expiry rather than removing it mechanically;
#   * exclude a very short expiry from the headline H estimate when
#     finite-difference robustness checks are unstable;
#   * use only headline local-quadratic skews with two-sided support.
H_SELECTION_MAX_TAU_DAYS <- 180L
H_SENSITIVITY_MAX_TAU_DAYS <- 365L

H_selection_audit <- atm_skew_tbl %>%
  left_join(
    atm_skew_cd_tbl %>%
      select(expiration, tau_days, tau, cd_skew_h0025, cd_skew_h005, cd_available_h0025, cd_available_h005, cd_sign_same_h0025, cd_sign_same_h005, cd_robust_for_headline, cd_robust_comment),
    by = c("expiration", "tau_days", "tau")
  ) %>%
  mutate(
    in_short_end_window = tau_days <= H_SELECTION_MAX_TAU_DAYS,
    passes_headline_estimator = headline_atm_skew,
    passes_cd_robustness = cd_robust_for_headline,
    used_for_headline_H = in_short_end_window & passes_headline_estimator & passes_cd_robustness & !very_short_maturity,
    used_for_sensitivity_including_1d = in_short_end_window & passes_headline_estimator & passes_cd_robustness,
    used_for_sensitivity_extended_1y = tau_days <= H_SENSITIVITY_MAX_TAU_DAYS & passes_headline_estimator & passes_cd_robustness & !very_short_maturity,
    H_selection_reason = case_when(
      used_for_headline_H ~ "used_for_headline_short_end_H",
      very_short_maturity & passes_headline_estimator & !passes_cd_robustness ~ "excluded_from_headline_H_after_CD_robustness_failure",
      !in_short_end_window & tau_days <= H_SENSITIVITY_MAX_TAU_DAYS & passes_headline_estimator & passes_cd_robustness ~ "outside_headline_short_end_window_used_only_in_extended_sensitivity",
      !passes_headline_estimator & good_atm_skew ~ "excluded_from_headline_H_lower_quality_skew_estimator",
      !good_atm_skew ~ "excluded_from_headline_H_failed_ATM_skew_quality",
      TRUE ~ "excluded_from_headline_H_other"
    )
  ) %>%
  arrange(tau_days)

headline_H_days <- H_selection_audit %>% filter(used_for_headline_H) %>% pull(tau_days)
calibration_good_expiries <- H_selection_audit %>%
  filter(good_atm_skew, tau_days <= H_SENSITIVITY_MAX_TAU_DAYS, !(very_short_maturity & !passes_cd_robustness)) %>%
  arrange(tau_days) %>%
  pull(tau_days)

if (FAST_TEST) {
  headline_H_days <- head(headline_H_days, 3)
  calibration_good_expiries <- head(calibration_good_expiries, 3)
}

loglog_headline <- fit_loglog_skew(
  atm_skew_tbl %>% filter(tau_days %in% headline_H_days),
  label = "headline_short_end_robust_for_calibration"
)

loglog_sens_include_1d_support <- fit_loglog_skew(
  atm_skew_tbl %>% filter(tau_days <= H_SELECTION_MAX_TAU_DAYS, headline_atm_skew),
  label = "sensitivity_short_end_support_only_including_1d"
)

loglog_sens_extended_1y <- fit_loglog_skew(
  atm_skew_tbl %>% filter(tau_days %in% (H_selection_audit %>% filter(used_for_sensitivity_extended_1y) %>% pull(tau_days))),
  label = "sensitivity_extended_to_1y_excluding_unstable_1d"
)

loglog_all_good_short_end <- fit_loglog_skew(
  atm_skew_tbl %>% filter(good_atm_skew, tau_days <= H_SELECTION_MAX_TAU_DAYS, tau_days > 2),
  label = "sensitivity_short_end_all_two_sided_excluding_1d"
)

atm_skew_loglog_summary <- bind_rows(
  loglog_headline$summary,
  loglog_sens_include_1d_support$summary,
  loglog_sens_extended_1y$summary,
  loglog_all_good_short_end$summary
)

atm_skew_loglog_points <- bind_rows(
  loglog_headline$points,
  loglog_sens_include_1d_support$points,
  loglog_sens_extended_1y$points,
  loglog_all_good_short_end$points
)

empirical_H_for_calibration <- atm_skew_loglog_summary %>%
  filter(estimator == "headline_short_end_robust_for_calibration") %>%
  mutate(
    selected_for_calibration = TRUE,
    selection_rule = paste0(
      "short-end maturities tau_days <= ", H_SELECTION_MAX_TAU_DAYS,
      "; local-quadratic ATM-skew estimator; two-sided near-ATM support; ",
      "central-difference signs consistent; one-day excluded after robustness failure if applicable"
    ),
    calibration_expiries_days = paste(calibration_good_expiries, collapse = ", "),
    interpretation = "empirical H used to fix alpha = H + 1/2 before rough Heston calibration"
  )

if (nrow(empirical_H_for_calibration) != 1L || !is.finite(empirical_H_for_calibration$H_for_calibration[1])) {
  stop("No valid headline empirical H estimate for calibration. Inspect phaseM2010_market_H_selection_audit.csv.")
}

# Calibration panel: broader than the skew-extraction panel, but restricted to
# explicitly audited maturities. The one-day expiry is not hard-excluded; in this
# data set it is excluded from the main calibration panel because its ATM-skew
# robustness check fails.
opts_panel <- opts_panel %>%
  mutate(
    good_expiry = tau_days %in% calibration_good_expiries,
    eligible_for_calibration = good_expiry &
      tau_days <= H_SENSITIVITY_MAX_TAU_DAYS &
      abs(k) <= 0.35 &
      rel_spread <= 1.0 &
      selected_is_otm &
      is.finite(iv_mkt),
    rank_abs_k = ave(abs(k), tau_days, FUN = function(x) rank(x, ties.method = "first")),
    holdout_flag = eligible_for_calibration & (rank_abs_k %% 5 == 0),
    calibration_flag = eligible_for_calibration & !holdout_flag,
    calibration_panel_label = ifelse(calibration_flag, "Yes", "No"),
    used_for_headline_H = tau_days %in% headline_H_days
  )

expiry_summary <- opts_panel %>%
  group_by(expiration, tau_days, tau) %>%
  summarise(
    listed_expiration = first(expiration),
    pricing_expiration = first(pricing_expiration),
    expiration_wday = first(expiration_wday),
    F = first(F),
    D = first(D),
    n_quotes_panel = n(),
    n_calibration = sum(calibration_flag),
    n_holdout = sum(holdout_flag),
    k_min = min(k, na.rm = TRUE),
    k_max = max(k, na.rm = TRUE),
    iv_min = min(iv_mkt, na.rm = TRUE),
    iv_max = max(iv_mkt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(atm_skew_tbl, by = c("expiration", "tau_days", "tau")) %>%
  left_join(
    H_selection_audit %>% select(expiration, tau_days, tau, used_for_headline_H, H_selection_reason, passes_cd_robustness),
    by = c("expiration", "tau_days", "tau")
  ) %>%
  arrange(tau_days)

audit_summary <- tibble(
  selected_date = target_date,
  S0 = S0,
  n_panel_quotes = nrow(opts_panel),
  n_calibration_quotes = sum(opts_panel$calibration_flag),
  n_holdout_quotes = sum(opts_panel$holdout_flag),
  n_expiries_with_good_atm_skew = sum(atm_skew_tbl$good_atm_skew, na.rm = TRUE),
  n_headline_H_maturities = length(headline_H_days),
  headline_H_maturities_days = paste(headline_H_days, collapse = ", "),
  calibration_expiries_days = paste(calibration_good_expiries, collapse = ", "),
  H_for_calibration = empirical_H_for_calibration$H_for_calibration[1],
  alpha_for_calibration = empirical_H_for_calibration$alpha_for_calibration[1],
  one_day_policy = "Very short maturities are audited, not mechanically excluded. The 1-day maturity is retained in the support audit but excluded from headline H and main calibration if its central-difference robustness check fails.",
  atm_skew_estimator = "headline local quadratic with two-sided near-ATM support; central differences are robustness diagnostics only",
  calibration_panel_policy = "broader OTM call-equivalent panel, restricted to explicitly audited maturities and abs(k) <= 0.35"
)

model_universe <- tibble(
  model = c("rough Heston", "classical Heston (alpha=1)", "Variance Gamma"),
  role = c("main", "benchmark", "benchmark"),
  mechanism = c("rough stochastic volatility", "Brownian stochastic volatility", "pure-jump exponential Levy"),
  pricing_route = c("fractional-Riccati CF + Lewis inversion", "alpha=1 Riccati/CF + Lewis inversion", "VG characteristic function + Lewis inversion")
)

write_csv(model_universe, file.path(out_dir, "phaseM2010_model_universe.csv"))
write_csv(opts_panel, file.path(out_dir, "phaseM2010_market_panel.csv"))
write_csv(opts_panel %>% filter(calibration_flag | holdout_flag),
          file.path(out_dir, "phaseM2010_market_calibration_panel.csv"))
write_csv(expiry_summary, file.path(out_dir, "phaseM2010_expiry_summary.csv"))
write_csv(atm_skew_tbl, file.path(out_dir, "phaseM2010_market_atm_skew.csv"))
write_csv(near_atm_support_audit, file.path(out_dir, "phaseM2010_market_atm_skew_extraction_audit.csv"))
write_csv(atm_skew_cd_tbl, file.path(out_dir, "phaseM2010_market_atm_skew_central_difference_check.csv"))
write_csv(H_selection_audit, file.path(out_dir, "phaseM2010_market_H_selection_audit.csv"))
write_csv(empirical_H_for_calibration, file.path(out_dir, "phaseM2010_empirical_H_for_calibration.csv"))
write_csv(atm_skew_loglog_summary, file.path(out_dir, "phaseM2010_market_atm_skew_loglog_summary.csv"))
write_csv(atm_skew_loglog_points, file.path(out_dir, "phaseM2010_market_atm_skew_loglog_points.csv"))
write_csv(audit_summary, file.path(out_dir, "phaseM2010_market_panel_audit_summary.csv"))

print(audit_summary)
print(expiry_summary %>% select(tau_days, n_quotes_panel, n_calibration, n_holdout, atm_iv, atm_skew, fit_type, good_atm_skew, headline_atm_skew, passes_cd_robustness, used_for_headline_H, H_selection_reason))
print(atm_skew_loglog_summary)
print(empirical_H_for_calibration)

# ============================================================
# 6. Plots
# ============================================================

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 8.5),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

plot_days <- head(calibration_good_expiries, ifelse(FAST_TEST, 3, 8))

p_smiles <- opts_panel %>%
  filter(tau_days %in% plot_days, abs(k) <= 0.25) %>%
  mutate(
    tau_label = factor(paste0(tau_days, "d"), levels = paste0(plot_days, "d")),
    calibration_panel_label = factor(calibration_panel_label, levels = c("No", "Yes"))
  ) %>%
  ggplot(aes(x = k, y = iv_mkt)) +
  geom_point(aes(shape = calibration_panel_label), size = 1.2, alpha = 0.85) +
  geom_line(aes(group = tau_label), linewidth = 0.35, alpha = 0.65) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey35") +
  facet_wrap(~ tau_label, scales = "free_y") +
  labs(
    x = "Forward log-moneyness, k = log(K/F)",
    y = "B&S implied volatility",
    shape = "Calibration panel"
  ) +
  plot_theme

plot_save(p_smiles, "plot_phaseM2010_market_iv_slices", width = 10, height = 6)

p_skew <- atm_skew_tbl %>%
  left_join(H_selection_audit %>% select(expiration, tau_days, tau, used_for_headline_H), by = c("expiration", "tau_days", "tau")) %>%
  filter(good_atm_skew, tau_days <= H_SENSITIVITY_MAX_TAU_DAYS) %>%
  mutate(headline_H_label = ifelse(used_for_headline_H, "Used for headline H", "Diagnostic only")) %>%
  ggplot(aes(x = tau_days, y = atm_skew)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(shape = headline_H_label), size = 2.0) +
  labs(
    x = "Maturity in calendar days",
    y = expression(paste("Signed ATM skew, ", Psi(tau))),
    shape = "H-estimation role"
  ) +
  plot_theme

plot_save(p_skew, "plot_phaseM2010_market_atm_skew", width = 7.2, height = 4.6)

headline_fit_points <- loglog_headline$points %>%
  arrange(tau_days) %>%
  mutate(tau_grid = tau, fitted_abs_atm_skew = exp(fitted_log_abs_atm_skew))

p_loglog <- atm_skew_tbl %>%
  left_join(H_selection_audit %>% select(expiration, tau_days, tau, used_for_headline_H, H_selection_reason), by = c("expiration", "tau_days", "tau")) %>%
  filter(good_atm_skew, tau_days <= H_SENSITIVITY_MAX_TAU_DAYS, is.finite(atm_skew), abs(atm_skew) > 0) %>%
  mutate(
    abs_atm_skew = abs(atm_skew),
    headline_H_label = ifelse(used_for_headline_H, "Used for headline H", "Diagnostic only")
  ) %>%
  ggplot(aes(x = tau, y = abs_atm_skew)) +
  geom_point(aes(shape = headline_H_label), size = 2.2) +
  geom_line(data = headline_fit_points, aes(x = tau_grid, y = fitted_abs_atm_skew), inherit.aes = FALSE, linewidth = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(paste("Time to maturity, ", tau, " (years, log scale)")),
    y = expression(paste("Absolute ATM skew, |", Psi(tau), "| (log scale)")),
    shape = "H-estimation role"
  ) +
  plot_theme

plot_save(p_loglog, "plot_phaseM2010_market_atm_skew_loglog", width = 7.2, height = 4.8)

figure_manifest <- tibble(
  file = paste0(c(
    "plot_phaseM2010_market_iv_slices",
    "plot_phaseM2010_market_atm_skew",
    "plot_phaseM2010_market_atm_skew_loglog"
  ), ".", PLOT_FORMAT),
  include_in_report = c(TRUE, TRUE, TRUE),
  description = c(
    "Clean 2010 SPX market B&S implied volatility slices used for empirical calibration.",
    "Market signed ATM-skew term structure with headline H-estimation maturities marked.",
    "Log-log market ATM-skew diagnostic used to estimate empirical H for rough Heston calibration."
  )
)
write_csv(figure_manifest, file.path(out_dir, "phaseM2010_figure_manifest.csv"))

cat("\n=== Phase M1 2010 market panel completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("Panel quotes:", nrow(opts_panel), "\n")
cat("Calibration quotes:", sum(opts_panel$calibration_flag), "\n")
cat("Holdout quotes:", sum(opts_panel$holdout_flag), "\n")
cat("Headline H maturities:", paste(headline_H_days, collapse = ", "), "\n")
cat("Calibration expiries:", paste(calibration_good_expiries, collapse = ", "), "\n")
cat("Empirical H for calibration:", empirical_H_for_calibration$H_for_calibration[1], "\n")
cat("Empirical alpha for calibration:", empirical_H_for_calibration$alpha_for_calibration[1], "\n")
cat("Next step: inspect Phase M1 outputs before running PhaseM_2010_rheston_calibration.R\n")
