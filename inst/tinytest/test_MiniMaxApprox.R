# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Test other components of minimaxApprox - Poly that have not been tested before
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
PP <- minimaxApprox(fn, -1, 1, 9, errType = "rel",
                    opts = list(maxiter = 100, convRatio = 1.5))
expect_false(PP$Warning)

# Should show at least one line of output due to show progress
expect_message(minimaxApprox(fn, -1, 1, 9,
                             opts = list(miniter = 0L, showProgress = TRUE)),
               "i: 1 E: ")

fn <- function(x) sin(x) + cos(x)
expect_warning(minimaxApprox(fn, -1, 1, 13),
               "All errors very near machine double precision.")
expect_warning(minimaxApprox(fn, -1, 1, 9, opts = list(maxiter = 0)),
               "Convergence not acheived in ")

PP <- suppressWarnings(minimaxApprox(fn, -1, 1, 13, opts = list(tol = 1e-14)))
expect_true(PP$Warning)

expect_warning(minimaxApprox(fn, -1, 1, 10, xi = 6),
               "Polynomial approximation uses Chebyeshev nodes for initial")

# Test other components of minimaxApprox - Rat that have not been tested before
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
dg <- c(2L, 1L)
expect_false(minimaxApprox(fn, -1, 1, dg)$Warning)

# Should show at least one line of output due to show progress
expect_message(minimaxApprox(fn, -1, 1, dg, opts = list(miniter = 2L,
                                                        showProgress = TRUE)),
               "i: 1 E: ")

dg <- c(3L, 3L)
# Test error message
expect_error(minimaxApprox(fn, -1, 1, dg, xi = chebNodes(5, -1, 1)),
             "Given the requested degrees for numerator and denominator")

expect_error(minimaxApprox(fn, -1, 1, 1:3),
             "Polynomial approximation takes one value for degree and")

# Test warning flag
expect_true(suppressWarnings(minimaxApprox(fn, -1, 1, c(5, 4),
                                           opts = list(maxiter = 6L))$Warning))

# Test informational message
fn <- function(x) exp(x) - 1
expect_message(minimaxApprox(fn, -0.15, 0.15, c(3, 4)),
               "All errors very near machine double precision. The solution")

# Test passing some parameters
expect_silent(minimaxApprox(fn, -0.15, 0.15, c(2, 2),
                            opts = list(convRatio = 1.05, tol = 1e-4)))

# Test denominator error message
expect_error(minimaxApprox(tan, 1, 2, c(2, 3)),
             "The 3 degree polynomial in the denominator has a zero at 1.57")

# Test print, plot, and coef methods
PP <- minimaxApprox(function(x) exp(x), 0, 1, 5, "abs")
expect_identical(unlist(coef(PP), use.names = FALSE), PP$a)
expect_stdout(print(PP))
expect_stdout(plot(PP))

RR <- minimaxApprox(function(x) exp(x), 0, 1, c(2L, 2L), "abs")
expect_identical(unlist(coef(RR)$a, use.names = FALSE), RR$a)
expect_identical(unlist(coef(RR)$b, use.names = FALSE), RR$b)
expect_stdout(print(RR))
expect_stdout(plot(RR))

RR <- suppressWarnings(minimaxApprox(exp, 0, 1, c(2L, 2L), "rel"))
expect_identical(unlist(coef(RR)$a, use.names = FALSE), RR$a)
expect_identical(unlist(coef(RR)$b, use.names = FALSE), RR$b)
expect_stdout(print(RR))
expect_stdout(plot(RR))
expect_stdout(plot(RR, ylim = c(-5e-6, 5e-6)))
