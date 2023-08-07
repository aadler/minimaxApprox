# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
x <- 1:3
a <- c(1, 1e-6, 1e-12)
expect_equal(minimaxApprox:::horner(x, a),
             minimaxApprox:::compensatedHorner(x, a),
             tolerance = tol)

expect_error(minimaxApprox:::hornerSum(x, a, c(a, 1)),
             "Error polynomials must be of same length.")
