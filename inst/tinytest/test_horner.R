# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7
x <- 1:3
a <- c(1, 1e-6, 1e-12)
aa <- matrix(rep(a, 2), ncol = 3)
testf <- function(x, a) sum(x ^ (0:2) * a)
control <- vapply(x, testf, double(1L), a = a)

## compensatedHorner
# Differences should be much less than tol - at the edges of machine precision
# or beyond.
expect_identical(minimaxApprox:::polyCalc(x, a), control)

## hornerSum
# Proper function
expect_silent(.Call(minimaxApprox:::hornerSum_c, as.double(x),
                    aa, NROW(aa), aa))

# Test error traps
expect_error(.Call(minimaxApprox:::hornerSum_c, as.double(x),
                   aa, NROW(aa), rbind(aa, aa)),
             "Error polynomials must be of same dimension.")

expect_error(.Call(minimaxApprox:::hornerSum_c, as.double(x[-1L]),
                   aa, NROW(aa), aa),
             "Polynomials must have same length as x.")

# Test 0 trap
aa <- double(0)
expect_identical(.Call(minimaxApprox:::hornerSum_c, as.double(x),
                       aa, NROW(aa), aa), rep(0, length(x)))
