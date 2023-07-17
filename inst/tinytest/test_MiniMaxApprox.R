# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- 1e-7

# Most tests of warnings and messages will perforce check internals too.

# Check Accuracy
## Rational. Based on Cody (1968) pp 250--251. Using weaker tolerance since
## taking values printed on paper.
controlA <- c(1.2655835, -0.65058499, 0.19786869)
controlB <- c(1, -0.064342748, -0.028851456)
controlX <- c(2, 2.0924, 2.3368, 2.6459, 2.9011, 3)
controlE <- 2.6934e-5
RR <- minimaxApprox(gamma, 2, 3, c(2L, 2L), TRUE, opts = list())
expect_equal(RR$a, controlA, tolerance = 5e-5)
expect_equal(RR$b, controlB, tolerance = 5e-5)
expect_equal(RR$x, controlX, tolerance = 5e-5)
expect_equal(RR$EE, controlE, tolerance = 5e-5)
expect_false(RR$Warning)

# Test warning flag in good case
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
## Polynomial
expect_false(minimaxApprox(fn, -1, 1, 9L)$Warning)
## Rational tested above

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
expect_error(minimaxApprox(sin, 0, pi / 4, c(1L, 1L), TRUE), errMess)

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
## Check error trap
errMess <- "This function only works with 'minimaxApprox' objects."
expect_error(minimaxEval(x, sin), errMess)
