test_that("mr_input validates dimensions", {
  expect_error(mr_input(1:3, rep(1, 3), matrix(1, 2, 2), matrix(1, 2, 2)))
  expect_error(mr_input(1:5, rep(-1, 5), matrix(1, 5, 2), matrix(1, 5, 2)))
})

test_that("simulated truth lies in the region", {
  sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 1)
  expect_true(searching_membership(sim, attr(sim, "true_beta"))$in_region)
  reg <- searching_region(sim)
  tb <- attr(sim, "true_beta")
  expect_true(all(tb >= reg$sections[, 1] & tb <= reg$sections[, 2]))
})

test_that("far-away points are excluded", {
  sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 1)
  expect_false(searching_membership(sim, c(20, 20))$in_region)
})

test_that("a minority-strong coordinate is unbounded", {
  sim <- simulate_mvmr(q = 2, p = 30, n_valid = 21, seed = 7)
  sim$se_X[6:30, 2] <- sim$se_X[6:30, 2] * 20
  b <- boundedness_rule(sim)
  expect_false(b$bounded[2])
  reg <- suppressWarnings(searching_region(sim))
  expect_true(any(is.infinite(reg$sections[2, ])))
})

test_that("boundedness threshold is floor(|S|/2)", {
  sim <- simulate_mvmr(q = 2, p = 31, n_valid = 22, seed = 3)
  S <- relevance_screen(sim)
  expect_equal(boundedness_rule(sim, S = S)$threshold[1], floor(length(S) / 2))
})

test_that("diagnostics and comparators run", {
  sim <- simulate_mvmr(q = 3, p = 60, n_valid = 40, seed = 2)
  d <- iv_diagnostics(sim)
  expect_s3_class(d, "iv_diagnostics")
  expect_equal(nrow(d$per_exposure), 3L)
  expect_true(is.numeric(mv_ivw(sim)$estimate))
  expect_true(is.numeric(tsht_valid_set(sim)$surplus))
})
