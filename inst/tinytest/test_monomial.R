# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# Test vanderMat
k <- 1:5
control <- matrix(c(rep(1, 5L), k, k ^ 2, k ^ 3, k ^ 4, k ^ 5, k ^ 6), ncol = 7)
expect_identical(minimaxApprox:::vanderMat(k, 6L), control)

# Test polyCalc
coeffs <-  c(2, 3.2, 4.6, -9.7, 0.1)
controlF <- function(x) {
  2 + 3.2 * x + 4.6 * x ^ 2 - 9.7 * x ^ 3 + 0.1 * x ^ 4
}

## Test scalar
x <- 3
expect_equal(minimaxApprox:::polyCalc(x, coeffs), controlF(x), tolerance = tol)
x <- 5
expect_equal(minimaxApprox:::polyCalc(x, coeffs), controlF(x), tolerance = tol)
x <- 1e-14
expect_equal(minimaxApprox:::polyCalc(x, coeffs), controlF(x), tolerance = tol)

## Test vectorized
expect_equal(minimaxApprox:::polyCalc(c(3, 5, 1e-14), coeffs),
             c(controlF(3), controlF(5), controlF(1e-14)), tolerance = tol)
