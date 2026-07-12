# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# chebMat
x <- c(-1.5, 0, 1.5)
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

# evalFunc using Chebyshev
## Uses a and x from immediately previous
## Polynomial
R <- list(a = a)
expect_equal(minimaxApprox:::evalFunc(x, R, "c"), control, tolerance = tol)

## Rational
b <- c(-0.5, -1, 2)
controlD <- c(-0.5 * 1 - (-0.5) + 2 * (2 * (-0.5) ^ 2 - 1),
              -0.5 * 1 - 1.5 + 2 * (2 * 1.5 ^ 2 - 1))
control <- control / controlD
R <- list(a = a, b = b)
expect_equal(minimaxApprox:::evalFunc(x, R, "c"), control, tolerance = tol)

# cheb2mon
## Polynomial
A <- minimaxApprox(exp, 0, 1, 4L)
B <- minimaxApprox(exp, 0, 1, 4L, basis = "m")
expect_equal(A$aMono, B$a, tolerance = tol)

## Rational
A <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L))
B <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L), basis = "m")

expect_equal(A$aMono, B$a, tolerance = tol)
expect_equal(A$bMono, B$b, tolerance = tol)

# --- F1 regression: chebCalc_c large-vector VLA/R_alloc fix -----------------
# Prior to the fix, src/Chebyshev.c's chebCalc_c allocated its working matrix
# as a C stack VLA (`double cMat[m * n]`); a several-million-point evaluation
# overflowed the C call stack and crashed the R session ("segfault from C
# stack overflow"). The fix (R_alloc) must handle this cleanly. One large call
# is used to keep runtime reasonable; a handful of spot-check expectations
# confirm correctness rather than just "didn't crash".
mmC <- minimaxApprox(exp, -1, 1, 5L, basis = "Chebyshev")
mmM <- minimaxApprox(exp, -1, 1, 5L, basis = "monomial")
xBig <- seq(0, 1, length.out = 5e6)
resC <- minimaxEval(xBig, mmC)
resM <- suppressMessages(minimaxEval(xBig, mmM))
expect_length(resC, 5e6)
expect_true(all(is.finite(resC)))
# Chebyshev and monomial representations of the same degree-5 approximation
# must agree to floating tolerance at every point on the grid.
expect_equal(resC, resM, tolerance = 1e-8)
# Spot-check endpoints and midpoint against direct exp() at loose tolerance
# (the polynomial is a degree-5 minimax fit, not exp itself).
expect_equal(resC[c(1, 2.5e6, 5e6)], exp(xBig[c(1, 2.5e6, 5e6)]),
             tolerance = 1e-3)

# --- F12 regression: 2^31-cell guard -----------------------------------
# chebMat_c and chebCalc_c share the same check_cell_cap() guard, so
# exercising it through chebMat_c (cheap: only the scalar degree argument
# need be large, no large vector allocation required) validates the shared
# code path used by both entry points.
expect_error(minimaxApprox:::chebMat(1:2, 2147483647),
             pattern = "exceeds the supported 2^31-cell limit",
             fixed = TRUE)
