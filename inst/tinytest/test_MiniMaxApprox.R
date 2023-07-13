# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Most tests of warnings and messages will perforce check internals too.

# Test warning flag in good case
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
## Polynomial
expect_false(minimaxApprox(fn, -1, 1, 9L)$Warning)
## Rational
expect_false(minimaxApprox(fn, -1, 1, c(2L, 1L))$Warning)

# Test showProgress. Also tests passing miniter
## Polynomial
opts <- list(miniter = 1L, showProgress = TRUE)
expect_message(minimaxApprox(fn, -1, 1, 9L, opts = opts), "i: 1 E: ")
## Rational
expect_message(minimaxApprox(fn, -1, 1, c(2L, 1L), opts = opts), "i: 1 E: ")

# Test passing some maxiter, convRatio, tol, and conviter. Also checks conviter
# overwrite.
fn <- function(x) exp(x) - 1
opts <- list(maxiter = 25L, convRatio = 1.01, tol = 1e-12, conviter = 50L)
## Polynomial
expect_silent(minimaxApprox(fn, -0.15, 0.15, 4L, opts = opts))
## Rational
expect_silent(minimaxApprox(fn, -0.15, 0.15, c(2L, 2L), opts = opts))

# Test maxiter warning and warning flag
opts <- list(maxiter = 2L)
wrnMess <- paste("Convergence to requested ratio and tolerance not acheived in",
                 "2 iterations.\nThe ratio is ")
## Polynomial
expect_warning(minimaxApprox(fn, -1, 1, 9L, opts = opts), wrnMess)
expect_true(suppressWarnings(minimaxApprox(fn, -1, 1, 9L, opts = opts)$Warning))
## Rational
dg <- c(2L, 2L)
expect_warning(minimaxApprox(fn, -1, 1, dg, opts = opts), wrnMess)
expect_true(suppressWarnings(minimaxApprox(fn, -1, 1, dg, opts = opts)$Warning))

# Test "very near machine double" warning message
wrnMess <- paste("All errors very near machine double precision. The solution",
                 "may not be optimal but should be best given the desired",
                 "precision and floating point limitations. Try a lower degree",
                 "if needed.")
## Polynomial
fn <- function(x) sin(x) + cos(x)
expect_warning(minimaxApprox(fn, -1, 1, 13L), wrnMess)
## Rational
fn <- function(x) exp(x) - 1
expect_warning(minimaxApprox(fn, -0.15, 0.15, c(3L, 4L)), wrnMess)

# Test consecutive unchanging check and message
fn <- function(x) exp(x) - 1
i <- 5L
opts <- list(conviter = i)
wrnMess <- paste(i, "succesive calculated errors were too close to each other",
                 "to warrant further iterations.\n")
## Polynomial
expect_warning(minimaxApprox(fn, -1, 1, 12L, opts = opts), wrnMess)
## Rational
expect_warning(minimaxApprox(fn, -1, 1, c(3L, 3L), opts = opts), wrnMess)

# Test function choosing basis x as 0 trap
errMess <- paste("Algorithm is choosing basis point where functional value is",
                 "0. Please approximate using absolute, and not relative",
                 "error.")
expect_error(minimaxApprox(sin, 0, pi / 4, c(1L, 1L), errType = 'rel'), errMess)

# Test passing incorrect degree (at minimaxApprox level)
errMess <- paste("Polynomial approximation takes one value for degree and",
                 "rational approximation takes a vector of two values for",
                 "numerator and denominator. Any other inputs are invalid.")
expect_error(minimaxApprox(fn, -1, 1, 1:3), errMess)

# Test passing xi
## Polynomial - Check that it is ignored
wrnMess <- paste("Polynomial approximation uses Chebyeshev nodes for initial",
                 "guess. Any passed xi is ignored.")
expect_warning(minimaxApprox(fn, -1, 1, 10L, xi = 6), wrnMess)
## Rational - Check that proper length is passed
errMess <- paste("Given the requested degrees for numerator and denominator,",
                 "the x-vector needs to have 8 elements.")
xi <- chebNodes(5L, -1, 1)
expect_error(minimaxApprox(fn, -1, 1, c(3L, 3L), xi = xi), errMess)

# Test checkDenom error message
expect_error(minimaxApprox(tan, 1, 2, c(2L, 3L)),
             "The 3 degree polynomial in the denominator has a zero at 1.57")

# Test evaluation function
x <- seq(0.1, 0.4, 0.025)
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$EE))
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L))
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$EE))

# Test print, plot, and coef methods
PP <- minimaxApprox(function(x) exp(x), 0, 1, 5L, "abs")
expect_identical(unlist(coef(PP), use.names = FALSE), PP$a)
expect_stdout(print(PP))
expect_stdout(plot(PP))

PP <- minimaxApprox(function(x) exp(x), 0, 1, 5L, "rel")
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
