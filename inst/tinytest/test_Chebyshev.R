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
x <- c(3, "A")
expect_error(minimaxApprox:::chebPoly(x, 3) , "is.numeric\\(x\\) is not TRUE")
expect_error(minimaxApprox:::chebPoly(3, -1), "k >= 0 is not TRUE")
