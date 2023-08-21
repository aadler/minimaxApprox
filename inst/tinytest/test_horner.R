# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
x <- 1:3
a <- c(1, 1e-6, 1e-12)
aa <- matrix(rep(a, 2), ncol = 3)
expect_equal(minimaxApprox:::hornerC(x, a),
             minimaxApprox:::compensatedHorner(x, a),
             tolerance = tol)

expect_error(minimaxApprox:::hornerSumC(x, a, c(a, 1)),
             "Error polynomials must be of same length.")

expect_silent(minimaxApprox:::hornerSumC(x, aa, aa))
expect_error(minimaxApprox:::hornerSumC(x[-1], aa, aa),
             "Polynomials must have same length as x.")
