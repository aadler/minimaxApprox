# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# chebMat
x <- c(-0.5, 0, 0.5)
nx <- length(x)
k <- 3
control <- matrix(c(rep(1, nx), x, 2 * x ^ 2 - 1, 4 * x ^ 3 - 3 * x), ncol = 4L,
                  byrow = FALSE)

expect_equal(minimaxApprox:::chebMat(x, k), control, tolerance = tol)

# chebCalc
a <- 2:5
x <- c(-0.5, 1.5)
control <- c(2 * 1 + 3 * -0.5 + 4 * (2 * (-0.5) ^ 2 - 1) +
               5 * (4 * (-0.5) ^ 3 - 3 * (-0.5)),
             2 * 1 + 3 * 1.5 + 4 * (2 * 1.5 ^ 2 - 1) +
               5 * (4 * 1.5 ^ 3 - 3 * 1.5))

expect_equal(minimaxApprox:::chebCalc(x, a), control, tolerance = tol)

# evalFuncCheb
## Uses a from immediately previous

## Polynomial
R <- list(a = a)
expect_equal(minimaxApprox:::evalFuncCheb(x, R), control, tolerance = tol)

## Rational
b <- c(-0.5, -1, 2)
controlD <- c(-0.5 * 1 - (-0.5) + 2 * (2 * (-0.5) ^ 2 - 1),
              -0.5 * 1 - 1.5 + 2 * (2 * 1.5 ^ 2 - 1))
control <- control / controlD
R <- list(a = a, b = b)

expect_equal(minimaxApprox:::evalFuncCheb(x, R), control, tolerance = tol)

# cheb2mon
## Polynomial
A <- minimaxApprox(exp, 0, 1, 4L)
B <- minimaxApprox(exp, 0, 1, 4L, basis = "c")

expect_equal(B$aMono, A$a, tolerance = tol)
## Rational
A <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L))
B <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L), basis = "c")

expect_equal(B$aMono, A$a, tolerance = tol)
expect_equal(B$bMono, A$b, tolerance = tol)
