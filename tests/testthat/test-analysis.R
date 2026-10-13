# Small synthetic fixtures only; no network, no dependence on the shipped CSVs
# (except the data-shape test, which is expected to fail, see below).

make_design <- function(n = 60, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  X <- cbind(1, x2, x1^3)
  list(X = X, x1 = x1, x2 = x2)
}

test_that("analyze_regression_model recovers known coefficients (noise-free)", {
  d <- make_design()
  theta <- c(2, -1.5, 0.5)
  y <- d$X %*% theta + rnorm(nrow(d$X), sd = 1e-8)
  fit <- analyze_regression_model(d$X, y, "toy")
  expect_equal(as.vector(fit$coefficients), theta, tolerance = 1e-5)
})

test_that("normal-equation OLS matches lm()", {
  d <- make_design()
  set.seed(2)
  y <- d$X %*% c(1, 2, 3) + rnorm(nrow(d$X))
  fit <- analyze_regression_model(d$X, y, "toy")
  ref <- coef(lm(y ~ d$X - 1))
  expect_equal(as.vector(fit$coefficients), unname(ref), tolerance = 1e-8)
})

test_that("RSS equals the sum of squared residuals", {
  d <- make_design()
  set.seed(3)
  y <- d$X %*% c(1, 2, 3) + rnorm(nrow(d$X))
  fit <- analyze_regression_model(d$X, y, "toy")
  expect_equal(fit$rss, sum(fit$residuals^2))
})

test_that("AIC and BIC follow their definitions and BIC penalises more for n > 7", {
  d <- make_design(n = 60)
  set.seed(4)
  y <- d$X %*% c(1, 2, 3) + rnorm(60)
  fit <- analyze_regression_model(d$X, y, "toy")
  k <- ncol(d$X); n <- 60
  expect_equal(fit$aic, 2 * k - 2 * fit$log_likelihood)
  expect_equal(fit$bic, k * log(n) - 2 * fit$log_likelihood)
  expect_gt(fit$bic - fit$aic, 0)   # k*(log(60) - 2) > 0
})

test_that("an irrelevant extra column lowers RSS but BIC prefers the smaller model", {
  d <- make_design(n = 80)
  set.seed(5)
  y <- d$X %*% c(1, 2, 3) + rnorm(80)
  small <- analyze_regression_model(d$X, y, "small")
  big <- analyze_regression_model(cbind(d$X, rnorm(80)), y, "big")
  expect_lte(big$rss, small$rss)          # extra column never raises RSS
  expect_lt(small$bic, big$bic)           # fixed seeds, so deterministic
})

test_that("a singular design matrix stops with a clear message", {
  d <- make_design()
  X <- cbind(d$X, d$X[, 2])               # duplicated column
  expect_error(analyze_regression_model(X, rnorm(nrow(X)), "dup"), "singular")
})

test_that("skewness_fn is 0 for symmetric data and signed for skewed data", {
  expect_equal(skewness_fn(c(-2, -1, 0, 1, 2)), 0)
  expect_gt(skewness_fn(c(1, 1, 1, 1, 10)), 0)
  expect_lt(skewness_fn(c(-10, 1, 1, 1, 1)), 0)
})

test_that("kurtosis_fn of a two-point symmetric sample is 1", {
  expect_equal(kurtosis_fn(c(-1, 1, -1, 1)), 1)
})

test_that("train/test: prediction on held-out rows has the right length and sane R^2", {
  d <- make_design(n = 100)
  set.seed(6)
  y <- as.vector(d$X %*% c(1, 2, 3) + rnorm(100, sd = 0.1))
  idx <- 1:floor(0.7 * 100)
  fit <- analyze_regression_model(d$X[idx, ], y[idx], "train")
  pred <- d$X[-idx, ] %*% fit$coefficients
  expect_length(pred, 30)
  r2 <- 1 - sum((y[-idx] - pred)^2) / sum((y[-idx] - mean(y[-idx]))^2)
  expect_gt(r2, 0.9)
})

test_that("calculate_ci returns one non-negative margin per row and widens with sigma^2", {
  d <- make_design(n = 50)
  m1 <- calculate_ci(d$X[1:10, ], 1)
  m2 <- calculate_ci(d$X[1:10, ], 4)
  expect_length(m1, 10)
  expect_true(all(m1 >= 0))
  expect_equal(m2, 2 * m1, tolerance = 1e-8)
})

test_that("ABC rejection logic: accepted draws concentrate near the true parameter", {
  # Same recipe as analysis.R section 6, on a toy problem with a vectorised loop.
  set.seed(7)
  n <- 100
  x <- rnorm(n); X <- cbind(1, x)
  y <- X %*% c(1, 2) + rnorm(n, sd = 0.5)
  fit <- analyze_regression_model(X, y, "toy")
  se <- sqrt(diag(sum(fit$residuals^2) / (n - 2) * solve(t(X) %*% X)))
  b <- fit$coefficients
  sims <- 20000
  p1 <- runif(sims, b[1] - 3 * se[1], b[1] + 3 * se[1])
  p2 <- runif(sims, b[2] - 3 * se[2], b[2] + 3 * se[2])
  rss_sim <- vapply(seq_len(sims), function(i) sum((y - X %*% c(p1[i], p2[i]))^2), 0)
  keep <- abs(rss_sim - fit$rss) < 0.5
  expect_gt(sum(keep), 0)
  expect_lt(abs(mean(p1[keep]) - b[1]), 3 * se[1])
  expect_lt(abs(mean(p2[keep]) - b[2]), 3 * se[2])
})

test_that("shipped CSVs have equal numbers of rows", {
  # EXPECTED FAILURE: time.csv=201, x.csv=198, y.csv=204 rows, so analysis.R's
  # data.frame() call fails. This is a real bug in the repo's data (review.md C1).
  root <- project_root()
  n <- sapply(c("time.csv", "x.csv", "y.csv"),
              function(f) nrow(read.csv(file.path(root, f), header = FALSE)))
  # expect_failure() passes while the bug exists and starts failing once the
  # data are fixed, which is the cue to turn this into a plain expect_equal().
  expect_failure(expect_equal(length(unique(n)), 1))
})

test_that("ABC accepted draws are narrower than the OLS standard error (under-dispersed posterior)", {
  # Characterisation of analysis.R section 6: tolerance 0.5 on |RSS_sim - RSS_obs|
  # accepts only a small ellipse around the OLS estimate, so the "posterior" spread
  # is smaller than the OLS standard error. Documented in review.md (A3).
  set.seed(8)
  n <- 100
  x <- rnorm(n); X <- cbind(1, x)
  y <- X %*% c(1, 2) + rnorm(n, sd = 1)
  fit <- analyze_regression_model(X, y, "toy")
  se <- sqrt(diag(sum(fit$residuals^2) / (n - 2) * solve(t(X) %*% X)))
  b <- fit$coefficients
  sims <- 20000
  p1 <- runif(sims, b[1] - 3 * se[1], b[1] + 3 * se[1])
  p2 <- runif(sims, b[2] - 3 * se[2], b[2] + 3 * se[2])
  rss_sim <- vapply(seq_len(sims), function(i) sum((y - X %*% c(p1[i], p2[i]))^2), 0)
  keep <- abs(rss_sim - fit$rss) < 0.5
  expect_gt(sum(keep), 50)
  expect_lt(sd(p2[keep]), se[2])
})

test_that("time labels in df_long stay aligned with their rows", {
  # EXPECTED FAILURE (bug A1 in review.md): analysis.R lines 85-88 build df_long with
  # pivot_longer() (row-interleaved) and then rep(df$time, 5) (block-repeated), so
  # the time column does not match the values. The two lines below are copied from
  # analysis.R because the script cannot be sourced without running everything.
  skip_if_not_installed("dplyr"); skip_if_not_installed("tidyr")
  library(dplyr); library(tidyr)
  df <- data.frame(time = c(0, 0.002, 0.004), x1 = c(1, 2, 3), y = c(10, 20, 30))
  df_long <- df %>% select(-time) %>%
    pivot_longer(everything(), names_to = "signal_type", values_to = "value") %>%
    mutate(time = rep(df$time, ncol(df) - 1))
  # x1 value 2 belongs to time 0.002; the script gives it time 0.004
  expect_failure(expect_equal(df_long$time[df_long$value == 2], 0.002))
})
