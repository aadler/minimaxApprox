# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convRatio = 1.000000001, tol = 1e-14)
# Test ratMat
x <- c(-0.4, 0.1, 0.3, 0.4)
E <- 0.5
y <- x ^ 2
altSign <- (-1) ^ (seq_along(x) - 1)
Evctr <- E * altSign
yvctr <- y + Evctr
control <- matrix(c(rep(1, 4L), x, x ^ 2, -(yvctr) * x, -(yvctr) * x ^ 2,
                    -altSign), nrow = length(x))
expect_identical(ratMat(x, E, y, 2L, 2L, TRUE), control)

# Test ratCoeffs
# If the function is a pure polynomial then coeffs should recover it exactly IF
# we pass a one-term zero-degree polynomial in the denominator!
# Therefore for, for 4 data points and the function x² + 2x + 3, the
# coefficients should return c(3, 2, 1) and the error should be 0.
fn <- function(x) x ^ 2 + 2 * x + 3
x <- seq(0, 2, length.out = 4)
control <- c(3, 2, 1)
RR <- ratCoeffs(x, 0, fn, 2L, 0L, TRUE)
expect_equal(RR$a, control, tolerance = tol)
expect_identical(RR$b, 1)
expect_equal(RR$E, 0, tolerance = tol)
