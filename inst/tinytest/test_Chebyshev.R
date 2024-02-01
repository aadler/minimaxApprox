# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# See https://oeis.org/wiki/Chebyshev_polynomials

x <- 3

# Test chebPoly
## Scalar
expect_equal(minimaxApprox:::chebPoly(x, 0), 1, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 1), x, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 2), 2 * x ^ 2 - 1, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 3), 4 * x ^ 3 - 3 * x, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 4), 8 * x ^ 4 - 8 * x ^ 2 + 1,
             tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 5), 16 * x ^ 5 - 20 * x ^ 3 + 5 * x,
             tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 6),
             32 * x ^ 6 - 48 * x ^ 4 + 18 * x ^ 2 - 1, tolerance = tol)

## Vector and in all three zones
x <- c(-10, -2, -0.5, -1e-10, 0, 1e-10, 0.5, 2, 10)
expect_equal(minimaxApprox:::chebPoly(x, 0), rep(1, length(x)), tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 1), x, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 2), 2 * x ^ 2 - 1, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 3), 4 * x ^ 3 - 3 * x, tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 4), 8 * x ^ 4 - 8 * x ^ 2 + 1,
             tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 5), 16 * x ^ 5 - 20 * x ^ 3 + 5 * x,
             tolerance = tol)
expect_equal(minimaxApprox:::chebPoly(x, 6),
             32 * x ^ 6 - 48 * x ^ 4 + 18 * x ^ 2 - 1, tolerance = tol)

## Error trapping
expect_error(minimaxApprox:::chebPoly("A", 3))
expect_error(minimaxApprox:::chebPoly(3, -1))
expect_error(minimaxApprox:::chebPoly(3, 1.1))

# chebMat
x <- c(-0.5, 0, 0.5)
k <- 3
control <- matrix(c(minimaxApprox:::chebPoly(x, 0L),
                    minimaxApprox:::chebPoly(x, 1L),
                    minimaxApprox:::chebPoly(x, 2L),
                    minimaxApprox:::chebPoly(x, 3L)),
                  ncol = 4L, byrow = FALSE)

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
