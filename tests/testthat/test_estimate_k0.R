################################################################################
#                     Tests of the file estimate_k0.R                          #
################################################################################

library(SmoothCurves)

test_that("estimate_k0_pilot is double", {
  expect_equal(estimate_k0_pilot(200), 21)
  expect_type(estimate_k0_pilot(200), "double")
})

test_that("estimate_k0_oracle is double", {
  expect_equal(estimate_k0_oracle(200, 0.5), 3)
  expect_type(estimate_k0_oracle(200, 0.5), "double")
})
