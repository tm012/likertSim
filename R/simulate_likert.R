# ---- Small-N safe Likert simulator ----
# Purpose: generate synthetic/sensitivity Likert datasets from a low-N real dataset.
# Important: synthetic rows should not be treated as real participants in inferential analyses.
# Assumptions:
# - FIRST column is the participant/coder ID.
# - All remaining columns are Likert items within K_min..K_max.

suppressPackageStartupMessages({
  has_psych  <- requireNamespace("psych",  quietly = TRUE)
  has_MASS   <- requireNamespace("MASS",   quietly = TRUE)
  has_Matrix <- requireNamespace("Matrix", quietly = TRUE)
  has_dplyr  <- requireNamespace("dplyr",  quietly = TRUE)
  has_tidyr  <- requireNamespace("tidyr",  quietly = TRUE)
})

# ---------- Utilities ----------
coerce_likert <- function(x, K_min, K_max){
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- trimws(x)
  x <- suppressWarnings(as.numeric(x))
  x <- as.integer(round(x))
  x[x < K_min] <- K_min
  x[x > K_max] <- K_max
  x
}

likert_midpoint <- function(K_min, K_max){
  as.integer(round((K_min + K_max) / 2))
}

build_thresholds_empirical <- function(x, K_min, K_max){
  x <- x[!is.na(x)]
  K_range <- K_max - K_min + 1L
  if (K_range < 2L) return(numeric(0))
  tab <- tabulate(x - K_min + 1L, nbins = K_range) + 0.5
  p <- tab / sum(tab)
  cum <- cumsum(p)
  eps <- 1e-6
  qnorm(pmin(pmax(cum[1:(K_range - 1L)], eps), 1 - eps))
}

build_thresholds_equal <- function(K_min, K_max){
  K_range <- K_max - K_min + 1L
  if (K_range < 2L) return(numeric(0))
  q <- seq(1 / K_range, (K_range - 1) / K_range, by = 1 / K_range)
  qnorm(q)
}

sanitize_corr <- function(R, J){
  if (J <= 0) return(matrix(numeric(0), 0, 0))
  if (J == 1) return(matrix(1, 1, 1))

  if (is.null(R) || length(R) == 0) R <- diag(J)
  R <- as.matrix(R)

  if (!all(dim(R) == c(J, J))) R <- diag(J)
  R[!is.finite(R)] <- 0
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R[R > 1] <- 1
  R[R < -1] <- -1
  R
}

make_pd <- function(R){
  R <- as.matrix(R)
  J <- ncol(R)
  R <- sanitize_corr(R, J)
  if (J <= 1) return(R)

  if (has_Matrix) {
    return(as.matrix(Matrix::nearPD(R, corr = TRUE)$mat))
  }

  eig <- eigen(R, symmetric = TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  R2 <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  d <- 1 / sqrt(diag(R2))
  D <- diag(d, nrow = length(d))
  R2 <- D %*% R2 %*% D
  diag(R2) <- 1
  R2
}

estimate_corr_poly_or_pearson <- function(items){
  J <- ncol(items)
  n_completeish <- sum(stats::complete.cases(items))
  if (J <= 1) {
    R <- matrix(1, 1, 1)
    attr(R, "method_used") <- "identity_one_item"
    return(R)
  }

  # Polychoric correlations are very unstable with very small N.
  # Try them only when there is a minimally useful number of complete rows.
  if (has_psych && n_completeish >= 10) {
    pc <- try(psych::polychoric(items, correct = TRUE, smooth = TRUE), silent = TRUE)
    if (!inherits(pc, "try-error") && is.list(pc) && !is.null(pc$rho)) {
      R <- sanitize_corr(pc$rho, J)
      attr(R, "method_used") <- "polychoric"
      return(R)
    }
  }

  R <- try(suppressWarnings(stats::cor(items, use = "pairwise.complete.obs")), silent = TRUE)
  if (inherits(R, "try-error")) R <- diag(J)
  R <- sanitize_corr(R, J)
  attr(R, "method_used") <- "pearson_or_identity_fallback"
  R
}

as_n_by_j <- function(x, N, J){
  if (J == 0) return(matrix(numeric(0), nrow = N, ncol = 0))
  if (is.null(dim(x))) return(matrix(x, nrow = N, ncol = J, byrow = FALSE))
  matrix(x, nrow = N, ncol = J)
}

# ---------- Fit model from real data ----------
fit_likert_model <- function(df, K_min, K_max,
                             threshold_strategy = c("empirical", "equal"),
                             keep_constant_items = TRUE){
  threshold_strategy <- match.arg(threshold_strategy)
  stopifnot(is.data.frame(df) || is.matrix(df))
  df <- as.data.frame(df)
  stopifnot(ncol(df) >= 2)
  stopifnot(!is.null(K_min), !is.null(K_max), K_max >= K_min)

  raw_item_names <- colnames(df)[-1]
  if (is.null(raw_item_names) || any(raw_item_names == "")) {
    raw_item_names <- paste0("item_", seq_len(ncol(df) - 1L))
  }

  items_all <- df[, -1, drop = FALSE]
  colnames(items_all) <- raw_item_names
  items_all[] <- lapply(items_all, coerce_likert, K_min = K_min, K_max = K_max)

  keep_variable <- vapply(items_all, function(x){
    ux <- unique(x[!is.na(x)])
    length(ux) > 1
  }, logical(1))

  constant_items <- names(items_all)[!keep_variable]
  constant_values <- integer(0)
  if (length(constant_items) > 0) {
    constant_values <- vapply(items_all[constant_items], function(x){
      vals <- x[!is.na(x)]
      if (length(vals) == 0) return(likert_midpoint(K_min, K_max))
      as.integer(round(stats::median(vals)))
    }, integer(1))
  }

  items <- items_all[, keep_variable, drop = FALSE]
  J <- ncol(items)

  if (J > 0) {
    cuts_list <- if (threshold_strategy == "empirical") {
      lapply(items, build_thresholds_empirical, K_min = K_min, K_max = K_max)
    } else {
      eq <- build_thresholds_equal(K_min, K_max)
      replicate(J, eq, simplify = FALSE)
    }

    R <- estimate_corr_poly_or_pearson(items)
    method_used <- attr(R, "method_used")
    R <- make_pd(R)

    mu_items <- sapply(items, function(x) mean(as.numeric(x), na.rm = TRUE))
    sd_items <- sapply(items, function(x) stats::sd(as.numeric(x), na.rm = TRUE))
    sd_items[!is.finite(sd_items) | sd_items == 0] <- 1e-6
  } else {
    cuts_list <- list()
    R <- matrix(numeric(0), 0, 0)
    method_used <- "no_variable_items"
    mu_items <- numeric(0)
    sd_items <- numeric(0)
  }

  list(
    K_min = K_min,
    K_max = K_max,
    R = R,
    cuts = cuts_list,
    variable_items = colnames(items),
    all_items = raw_item_names,
    corr_method = method_used,
    mu = mu_items,
    sd = sd_items,
    constant_items = constant_items,
    constant_values = constant_values,
    keep_constant_items = keep_constant_items
  )
}

# ---------- Simulate rows ----------
simulate_likert <- function(model, N_sim,
                            id_prefix = "syn",
                            noise_factor = 1.0,
                            bias = 0.0,
                            corr_alpha = 0.5,
                            seed = NA,
                            mode = c("meansd", "threshold")){
  mode <- match.arg(mode)
  stopifnot(is.list(model), all(c("K_min", "K_max", "R", "variable_items", "all_items") %in% names(model)))
  stopifnot(length(N_sim) == 1, is.finite(N_sim), N_sim > 0)
  N_sim <- as.integer(N_sim)
  if (!is.null(seed) && length(seed) == 1 && !is.na(seed)) set.seed(seed)

  J <- length(model$variable_items)
  X_var <- matrix(numeric(0), nrow = N_sim, ncol = 0)

  if (J > 0) {
    a <- max(0, min(1, corr_alpha))
    R_stressed <- (1 - a) * model$R + a * diag(J)
    R_stressed <- make_pd(R_stressed)

    if (mode == "threshold") {
      Sigma <- R_stressed * (abs(noise_factor)^2)
      bias_vec <- rep(bias, J)

      if (has_MASS) {
        Z <- MASS::mvrnorm(n = N_sim, mu = bias_vec, Sigma = Sigma, empirical = FALSE)
        Z <- as_n_by_j(Z, N_sim, J)
      } else {
        eig <- eigen(Sigma, symmetric = TRUE)
        L <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8)), nrow = length(eig$values)) %*% t(eig$vectors)
        Z <- matrix(rnorm(N_sim * J), N_sim, J) %*% L
        Z <- sweep(Z, 2, bias_vec, FUN = "+")
      }

      X_var <- matrix(NA_integer_, nrow = N_sim, ncol = J)
      for (j in seq_len(J)) {
        X_var[, j] <- model$K_min + findInterval(Z[, j], model$cuts[[j]])
      }
    } else {
      mu <- as.numeric(model$mu)
      sd <- as.numeric(model$sd)
      mu_use <- mu + bias
      Dsd <- diag(pmax(sd, 1e-8), nrow = J, ncol = J)
      Sigma <- Dsd %*% R_stressed %*% Dsd
      Sigma <- Sigma * (abs(noise_factor)^2)

      if (has_MASS) {
        Yc <- MASS::mvrnorm(n = N_sim, mu = mu_use, Sigma = Sigma, empirical = FALSE)
        Yc <- as_n_by_j(Yc, N_sim, J)
      } else {
        eig <- eigen(Sigma, symmetric = TRUE)
        L <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8)), nrow = length(eig$values)) %*% t(eig$vectors)
        Yc <- matrix(rnorm(N_sim * J), N_sim, J) %*% L
        Yc <- sweep(Yc, 2, mu_use, FUN = "+")
      }

      X_var <- round(Yc)
      X_var[X_var < model$K_min] <- model$K_min
      X_var[X_var > model$K_max] <- model$K_max
    }
    colnames(X_var) <- model$variable_items
  }

  out_items <- as.data.frame(matrix(NA_integer_, nrow = N_sim, ncol = length(model$all_items)))
  colnames(out_items) <- model$all_items

  if (J > 0) {
    out_items[model$variable_items] <- as.data.frame(X_var, check.names = FALSE)
  }

  if (length(model$constant_items) > 0 && isTRUE(model$keep_constant_items)) {
    for (nm in model$constant_items) {
      out_items[[nm]] <- rep(model$constant_values[[nm]], N_sim)
    }
  }

  data.frame(coder_id = sprintf("%s_%05d", id_prefix, seq_len(N_sim)), out_items, check.names = FALSE)
}

mean_sd_by_item <- function(df, K_min = NULL, K_max = NULL){
  items <- df[, -1, drop = FALSE]
  if (!is.null(K_min) && !is.null(K_max)) {
    items[] <- lapply(items, coerce_likert, K_min = K_min, K_max = K_max)
  }
  m <- sapply(items, function(x) mean(as.numeric(x), na.rm = TRUE))
  s <- sapply(items, function(x) stats::sd(as.numeric(x), na.rm = TRUE))
  data.frame(item = names(m), mean = as.numeric(m), sd = as.numeric(s), row.names = NULL)
}

simulate_and_compare <- function(df_real,
                                 K_min,
                                 K_max,
                                 N_sim = nrow(df_real),
                                 noise_factor = 1.15,
                                 bias = 0.0,
                                 corr_alpha = 0.5,
                                 seed_default = 123,
                                 seed_stress = 456,
                                 id_prefix_default = "syn_base",
                                 id_prefix_stress = "syn_stress",
                                 mode = c("meansd", "threshold"),
                                 auto_small_switch = TRUE,
                                 small_n_cutoff = 15,
                                 threshold_strategy = c("empirical", "equal")){
  mode <- match.arg(mode)
  threshold_strategy <- match.arg(threshold_strategy)
  df_real <- as.data.frame(df_real)
  stopifnot(ncol(df_real) >= 2)

  n_raters <- nrow(df_real)
  chosen_mode <- mode
  if (auto_small_switch && n_raters < small_n_cutoff) {
    chosen_mode <- "meansd"
    message(sprintf("Small N detected (N=%d < %d). Using mean-SD mode with correlation shrinkage.",
                    n_raters, small_n_cutoff))
  }

  model <- fit_likert_model(df_real, K_min = K_min, K_max = K_max,
                            threshold_strategy = threshold_strategy)

  sim_default <- simulate_likert(model, N_sim = N_sim,
                                 id_prefix = id_prefix_default,
                                 noise_factor = 1.0,
                                 bias = 0.0,
                                 corr_alpha = corr_alpha,
                                 seed = seed_default,
                                 mode = chosen_mode)

  sim_stress <- simulate_likert(model, N_sim = N_sim,
                                id_prefix = id_prefix_stress,
                                noise_factor = noise_factor,
                                bias = bias,
                                corr_alpha = corr_alpha,
                                seed = seed_stress,
                                mode = chosen_mode)

  real_stats <- mean_sd_by_item(df_real, K_min = K_min, K_max = K_max); real_stats$dataset <- "real"
  def_stats  <- mean_sd_by_item(sim_default, K_min = K_min, K_max = K_max); def_stats$dataset <- "synthetic_base"
  str_stats  <- mean_sd_by_item(sim_stress, K_min = K_min, K_max = K_max); str_stats$dataset <- "synthetic_stress"

  comp <- rbind(real_stats, def_stats, str_stats)
  comp <- comp[, c("item", "dataset", "mean", "sd")]

  if (has_dplyr && has_tidyr) {
    comp <- comp |>
      dplyr::mutate(dataset = factor(dataset, levels = c("real", "synthetic_base", "synthetic_stress"))) |>
      dplyr::arrange(item, dataset)
  }

  sat_msg <- character(0)
  ds_list <- list(synthetic_base = sim_default, synthetic_stress = sim_stress)
  for (nm in names(ds_list)) {
    cur <- ds_list[[nm]]
    items_only <- cur[, -1, drop = FALSE]
    top_rate <- sapply(items_only, function(x) mean(as.integer(x) == K_max, na.rm = TRUE))
    bot_rate <- sapply(items_only, function(x) mean(as.integer(x) == K_min, na.rm = TRUE))
    sat_idx <- which(top_rate > 0.98 | bot_rate > 0.98)
    if (length(sat_idx) > 0) {
      sat_items <- names(items_only)[sat_idx]
      sat_msg <- c(sat_msg, sprintf("[%s] saturation on: %s", nm, paste(sat_items, collapse = ", ")))
    }
  }

  list(
    model = model,
    synthetic_base = sim_default,
    synthetic_stress = sim_stress,
    comparison = comp,
    warnings = sat_msg,
    mode_used = chosen_mode,
    threshold_strategy = threshold_strategy,
    note = "Use synthetic data for sensitivity/pipeline checks only. Do not treat synthetic rows as real participants."
  )
}

# Convenience wrapper for low-N projects.
make_synthetic_smallN <- function(df_real,
                                  K_min,
                                  K_max,
                                  target_N = 100,
                                  seed = 123,
                                  corr_alpha = 0.5,
                                  noise_factor = 1.15,
                                  bias = 0.0){
  simulate_and_compare(
    df_real = df_real,
    K_min = K_min,
    K_max = K_max,
    N_sim = target_N,
    mode = "meansd",
    auto_small_switch = TRUE,
    small_n_cutoff = 15,
    corr_alpha = corr_alpha,
    noise_factor = noise_factor,
    bias = bias,
    seed_default = seed,
    seed_stress = seed + 1,
    threshold_strategy = "equal"
  )
}
