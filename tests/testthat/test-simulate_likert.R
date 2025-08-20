test_that("fit/simulate basic flow", {
  set.seed(1)
  df <- data.frame(
    coder_id = paste0("id_", 1:20),
    A = sample(1:5, 20, TRUE),
    B = sample(1:5, 20, TRUE),
    C = sample(1:5, 20, TRUE)
  )
  model <- fit_likert_model(df, K_min = 1, K_max = 5)
  sim   <- simulate_likert(model, N_sim = 50, mode = "threshold", seed = 123)
  expect_equal(nrow(sim), 50)
  expect_equal(ncol(sim), 1 + length(model$items))
})
