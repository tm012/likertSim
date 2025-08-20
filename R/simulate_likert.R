#' Coerce a vector into a bounded Likert range
#' @keywords internal
coerce_likert <- function(x, K_min, K_max){
  x <- suppressWarnings(as.integer(round(x)))
  x[x < K_min] <- K_min
  x[x > K_max] <- K_max
  x
}

#' Build empirical thresholds for ordinal cutting (Laplace-smoothed)
#' @keywords internal
build_thresholds_empirical <- function(x, K_min, K_max){
  x <- x[!is.na(x)]
  K_range <- K_max - K_min + 1L
  if (K_range < 2L) return(numeric(0))
  tab <- tabulate(x - K_min + 1L, nbins = K_range) + 0.5
  p   <- tab / sum(tab)
  cum <- cumsum(p)
  eps <- 1e-6
  qnorm(pmin(pmax(cum[1:(K_range - 1L)], eps), 1 - eps))
}

#' Build equal-probability thresholds on N(0,1)
#' @keywords internal
build_thresholds_equal <- function(K_min, K_max){
  K_range <- K_max - K_min + 1L
  if (K_range < 2L) return(numeric(0))
  q <- seq(1 / K_range, (K_range - 1) / K_range, by = 1 / K_range)
  qnorm(q)
}

#' Estimate item correlation with polychoric (fallback to Pearson)
#' @keywords internal
estimate_corr_poly_or_pearson <- function(items){
  if (requireNamespace("psych", quietly = TRUE)) {
    pc <- try(psych::polychoric(items, correct = TRUE, smooth = TRUE), silent = TRUE)
    if (!inherits(pc, "try-error") && is.list(pc) && !is.null(pc$rho)) {
      attr(pc$rho, "method_used") <- "polychoric"
      return(pc$rho)
    }
  }
  R <- suppressWarnings(stats::cor(items, use = "pairwise.complete.obs"))
  attr(R, "method_used") <- "pearson_fallback"
  R
}

#' Make a correlation matrix positive-definite
#' @keywords internal
make_pd <- function(R){
  if (requireNamespace("Matrix", quietly = TRUE)) {
    return(as.matrix(Matrix::nearPD(R, corr = TRUE)$mat))
  }
  eig <- eigen(R, symmetric = TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  R2 <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  d  <- 1 / sqrt(diag(R2)); D <- diag(d)
  D %*% R2 %*% D
}

#' Fit a Likert simulation model from real data
#'
#' @param df data.frame, first column is `coder_id`, remaining columns are Likert items.
#' @param K_min integer, minimum Likert category (inclusive).
#' @param K_max integer, maximum Likert category (inclusive).
#' @param threshold_strategy "empirical" (learn per-item cuts) or "equal".
#' @return A list model with K_min, K_max, item names, thresholds, correlation,
#'   means/SDs, and diagnostics.
#' @export
fit_likert_model <- function(df, K_min, K_max,
                             threshold_strategy = c("empirical","equal")){
  threshold_strategy <- match.arg(threshold_strategy)
  stopifnot(is.data.frame(df) || is.matrix(df))
  df <- as.data.frame(df)
  stopifnot(ncol(df) >= 2, !is.null(K_min), !is.null(K_max), K_max >= K_min)

  items <- df[, -1, drop = FALSE]
  items[] <- lapply(items, coerce_likert, K_min = K_min, K_max = K_max)

  keep <- vapply(items, function(x){
    ux <- unique(x[!is.na(x)])
    length(ux) > 1
  }, logical(1))
  dropped_items <- names(items)[!keep]
  if (any(!keep)) {
    warning(sprintf("Dropped %d item(s) with no variance after coercion: %s",
                    sum(!keep), paste(dropped_items, collapse = ", ")))
    items <- items[, keep, drop = FALSE]
  }
  if (ncol(items) < 1) stop("No usable items remain after dropping zero-variance columns.")

  cuts_list <- if (threshold_strategy == "empirical") {
    lapply(items, build_thresholds_empirical, K_min = K_min, K_max = K_max)
  } else {
    eq <- build_thresholds_equal(K_min, K_max)
    replicate(ncol(items), eq, simplify = FALSE)
  }

  R <- estimate_corr_poly_or_pearson(items)
  method_used <- attr(R, "method_used")
  R <- make_pd(R)

  mu_items <- sapply(items, function(x) mean(as.numeric(x), na.rm = TRUE))
  sd_items <- sapply(items, function(x) stats::sd(as.numeric(x),   na.rm = TRUE))
  sd_items[!is.finite(sd_items) | sd_items == 0] <- 1e-6

  list(
    K_min = K_min,
    K_max = K_max,
    R     = R,
    cuts  = cuts_list,
    items = colnames(items),
    corr_method = method_used,
    mu    = mu_items,
    sd    = sd_items,
    dropped_items = dropped_items
  )
}

#' Simulate Likert responses from a fitted model
#'
#' @param model A model from \code{fit_likert_model()}.
#' @param N_sim Number of rows to simulate.
#' @param id_prefix Prefix for synthetic coder IDs.
#' @param noise_factor Scalar SD multiplier (dispersion stress).
#' @param bias Scalar mean shift (leniency/severity).
#' @param corr_alpha Blend toward independence: 0 keeps R, 1 = identity.
#' @param seed Optional RNG seed.
#' @param mode "threshold" (ordinal latent + cuts) or "meansd" (match means/SDs).
#' @return A data.frame with `coder_id` + item columns.
#' @export
simulate_likert <- function(model, N_sim,
                            id_prefix    = "sim",
                            noise_factor = 1.0,
                            bias         = 0.0,
                            corr_alpha   = 0.0,
                            seed         = NA,
                            mode         = c("threshold","meansd")) {
  mode <- match.arg(mode)
  stopifnot(is.list(model), all(c("K_min","K_max","R","items") %in% names(model)))
  if (!is.na(seed)) set.seed(seed)
  J <- length(model$items)

  a <- max(0, min(1, corr_alpha))
  R_stressed <- (1 - a) * model$R + a * diag(J)
  R_stressed <- make_pd(R_stressed)

  if (mode == "threshold") {
    Sigma <- R_stressed * (noise_factor^2)
    bias_vec <- rep(bias, J)

    if (requireNamespace("MASS", quietly = TRUE)) {
      Z <- MASS::mvrnorm(n = N_sim, mu = bias_vec, Sigma = Sigma, empirical = FALSE)
    } else {
      eig <- eigen(Sigma, symmetric = TRUE)
      L <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8))) %*% t(eig$vectors)
      Z <- matrix(rnorm(N_sim * J), N_sim, J) %*% L
      Z <- sweep(Z, 2, bias_vec, FUN = "+")
    }

    X <- matrix(NA_integer_, nrow = N_sim, ncol = J)
    for (j in seq_len(J)) {
      X[, j] <- model$K_min + findInterval(Z[, j], model$cuts[[j]])
    }

  } else {
    mu <- as.numeric(model$mu)
    sd <- as.numeric(model$sd)

    mu_use <- mu + bias
    Dsd <- diag(pmax(sd, 1e-8))
    Sigma <- Dsd %*% R_stressed %*% Dsd
    Sigma <- Sigma * (noise_factor^2)

    if (requireNamespace("MASS", quietly = TRUE)) {
      Yc <- MASS::mvrnorm(n = N_sim, mu = mu_use, Sigma = Sigma, empirical = FALSE)
    } else {
      eig <- eigen(Sigma, symmetric = TRUE)
      L <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8))) %*% t(eig$vectors)
      Yc <- matrix(rnorm(N_sim * J), N_sim, J) %*% L
      Yc <- sweep(Yc, 2, mu_use, FUN = "+")
    }

    X <- round(Yc)
    X[X < model$K_min] <- model$K_min
    X[X > model$K_max] <- model$K_max
  }

  out <- as.data.frame(X)
  colnames(out) <- model$items
  data.frame(coder_id = sprintf("%s_%05d", id_prefix, seq_len(N_sim)), out, check.names = FALSE)
}

#' Per-item mean/SD summary (first column is ID)
#' @keywords internal
mean_sd_by_item <- function(df){
  items <- df[, -1, drop = FALSE]
  m <- sapply(items, function(x) mean(as.numeric(x), na.rm = TRUE))
  s <- sapply(items, function(x) stats::sd(as.numeric(x),   na.rm = TRUE))
  data.frame(item = names(m), mean = as.numeric(m), sd = as.numeric(s), row.names = NULL)
}

#' Fit, simulate (default vs stressed), and compare summaries
#'
#' @param df_real data.frame; first col = coder_id, others = Likert items.
#' @param K_min,K_max scale bounds.
#' @param N_sim rows to simulate in each dataset.
#' @param noise_factor,bias,corr_alpha stress knobs for `sim_stress`.
#' @param mode "threshold" or "meansd".
#' @param auto_small_switch switch to "meansd" if raters < `small_n_cutoff`.
#' @param small_n_cutoff threshold for auto-switch.
#' @param threshold_strategy "empirical" or "equal" for thresholds.
#' @return A list with model, sim_default, sim_stress, a tidy comparison, warnings, and metadata.
#' @export
simulate_and_compare <- function(df_real,
                                 K_min,
                                 K_max,
                                 N_sim         = nrow(df_real),
                                 noise_factor  = 1.0,
                                 bias          = 0.0,
                                 corr_alpha    = 0.0,
                                 seed_default  = 123,
                                 seed_stress   = 456,
                                 id_prefix_default = "sim",
                                 id_prefix_stress  = "stress",
                                 mode          = c("threshold","meansd"),
                                 auto_small_switch = TRUE,
                                 small_n_cutoff    = 15,
                                 threshold_strategy = c("empirical","equal")) {

  mode <- match.arg(mode)
  threshold_strategy <- match.arg(threshold_strategy)

  df_real <- as.data.frame(df_real)
  stopifnot(ncol(df_real) >= 2)

  model <- fit_likert_model(df_real, K_min = K_min, K_max = K_max,
                            threshold_strategy = threshold_strategy)

  n_raters <- nrow(df_real)
  chosen_mode <- mode
  if (auto_small_switch && identical(mode, "threshold") && n_raters < small_n_cutoff) {
    message(sprintf("Few raters detected (N=%d < %d). Switching to mean–SD mode.",
                    n_raters, small_n_cutoff))
    chosen_mode <- "meansd"
  }

  sim_default <- simulate_likert(model, N_sim = N_sim,
                                 id_prefix   = id_prefix_default,
                                 noise_factor = 1.0,
                                 bias = 0.0,
                                 corr_alpha = 0.0,
                                 seed        = seed_default,
                                 mode        = chosen_mode)

  sim_stress  <- simulate_likert(model, N_sim = N_sim,
                                 id_prefix   = id_prefix_stress,
                                 noise_factor = noise_factor,
                                 bias = bias,
                                 corr_alpha = corr_alpha,
                                 seed        = seed_stress,
                                 mode        = chosen_mode)

  real_stats <- mean_sd_by_item(df_real);      real_stats$dataset <- "real"
  def_stats  <- mean_sd_by_item(sim_default);  def_stats$dataset  <- "sim_default"
  str_stats  <- mean_sd_by_item(sim_stress);   str_stats$dataset  <- "sim_stress"

  comp <- rbind(real_stats, def_stats, str_stats)
  comp <- comp[, c("item","dataset","mean","sd")]
  if (requireNamespace("dplyr", quietly = TRUE) &&
      requireNamespace("tidyr", quietly = TRUE)) {
    comp <- comp |>
      dplyr::mutate(dataset = factor(dataset, levels = c("real","sim_default","sim_stress"))) |>
      dplyr::arrange(item, dataset)
  }

  sat_msg <- character(0)
  ds_list <- list(sim_default = sim_default, sim_stress = sim_stress)
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
    sim_default = sim_default,
    sim_stress  = sim_stress,
    comparison  = comp,
    warnings    = sat_msg,
    mode_used   = chosen_mode,
    threshold_strategy = threshold_strategy
  )
}
