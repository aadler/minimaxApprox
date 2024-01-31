# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# Most tests of warnings and messages will perforce check internals too.

# Check Accuracy and lack of warning flag when converged
## Rational 1: Based on Fraser & Hart (1962) p. 403 Table 2
controlA <- c(0.99999998510030375, 0.601781180619504719,
              0.186144903531821877, 0.0687440518995425058)
controlB <- c(1, 1.17899457599300466, -0.122321311431112167,
              -0.260995866188425578, 0.0609927504415305534)
controlE <- 1e-6
fn <- function(x) gamma(x + 1)
RR <- minimaxApprox(fn, 0, 1, c(3L, 4L), relErr = FALSE)
expect_equal(RR$a, controlA, tolerance = tol)
expect_equal(RR$b, controlB, tolerance = tol)
expect_true(RR$ExpErr <= controlE)
expect_false(RR$Warning)
## Rational 2: Based on Cody (1968) pp 250--251. Using weaker tolerance since
## taking values printed on paper.
controlA <- c(1.2655835, -0.65058499, 0.19786869)
controlB <- c(1, -0.064342748, -0.028851456)
controlX <- c(2, 2.0924, 2.3368, 2.6459, 2.9011, 3)
controlE <- 2.6934e-5
RR <- minimaxApprox(gamma, 2, 3, c(2L, 2L), relErr = TRUE, opts = list())
expect_equal(RR$a, controlA, tolerance = 5e-6)
expect_equal(RR$b, controlB, tolerance = 5e-6)
expect_equivalent(RR$Extrema, controlX, tolerance = 5e-5)
expect_equal(RR$ExpErr, controlE, tolerance = 5e-5)
expect_false(RR$Warning)

## Rational 3: Based on DLMF 3.11.19 https://dlmf.nist.gov/3.11#iii
# Difference on Windows machine is roughly 8.13e-6
controlA <- c(0.99999998917854, -0.34038938209347, -0.18915483763222,
              0.06658319420166)
controlB <- c(1, -0.34039052338838, 0.06086501629812, -0.01864476809090)
fn <- function(x) besselJ(x, nu = 0)
b0 <- 0.893576966279167522
RR <- minimaxApprox(fn, 0, b0, c(3L, 3L))
expect_equal(RR$a, controlA, tolerance = 1e-5)
expect_equal(RR$b, controlB, tolerance = 1e-5)
expect_false(RR$Warning)

# Test incorrect basis for analysis
errMsg <- "Must select either 'Cheb'yshev or 'mono'mial basis for analysis."
expect_error(minimaxApprox(exp, 0, 1, 0, basis = "x"), errMsg)
expect_error(minimaxApprox(exp, 0, 1, 0, basis = 4), errMsg)

# Test passing length 0 for polynomials or rationals
expect_silent(minimaxApprox(exp, 0, 1, 0))
expect_silent(minimaxApprox(exp, 0, 1, c(0, 1)))
expect_silent(minimaxApprox(exp, 0, 1, c(1, 0)))
expect_silent(minimaxApprox(exp, 0, 1, c(0, 0)))
expect_identical(minimaxApprox(exp, 0, 1, c(0, 0))$a,
                 minimaxApprox(exp, 0, 1, 0)$a)
expect_identical(minimaxApprox(exp, 0, 1, c(0, 0))$b, 1)
expect_identical(minimaxApprox(exp, 0, 1, c(3, 0))$a,
                 minimaxApprox(exp, 0, 1, 3)$a)
expect_identical(minimaxApprox(exp, 0, 1, c(3, 0))$b, 1)

# Test negative and integer trap
errMsg <- "Polynomial degrees must be integers of least 0 (constant)."
## Polynomial
expect_error(minimaxApprox(exp, 0, 1, 0.2), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, -1L), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, -2), errMsg, fixed = TRUE)
## Rational
expect_error(minimaxApprox(exp, 0, 1, c(1, -2)), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, c(1.2, 2)), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, c(-1.2, 8.01)), errMsg, fixed = TRUE)

# Test trap for relErr
errMsg <- paste("Relative Error must be a logical value. Default FALSE",
                 "returns absolute error.")
## Polynomial
expect_error(minimaxApprox(exp, -1, 1, 9L, "abs"), errMsg)
## Rational
expect_error(minimaxApprox(exp, -1, 1, c(3L, 3L), "abs"), errMsg)

# Test showProgress; also tests passing miniter.
## Polynomial
fn <- function(x) exp(x) - 1
opts <- list(miniter = 1L, showProgress = TRUE)
expect_message(minimaxApprox(fn, -1, 1, 9L, opts = opts), "i: 1 E: ")
## Rational
expect_message(minimaxApprox(fn, -1, 1, c(2L, 1L), opts = opts), "i: 1 E: ")

# Test passing some maxiter, convrat, tol, and conviter. Also checks conviter
# overwrite.
fn <- function(x) exp(x) - 1
opts <- list(maxiter = 25L, convrat = 1.01, tol = 1e-12, conviter = 50L)
## Polynomial
expect_silent(minimaxApprox(fn, -0.15, 0.15, 4L, opts = opts))
## Rational
expect_silent(minimaxApprox(fn, -0.15, 0.15, c(2L, 2L), opts = opts))

# Test maxiter warning and warning flag
opts <- list(maxiter = 2L)
wrnMess <- paste("Convergence to requested ratio and tolerance not achieved in",
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
                 "may not be optimal given floating point limitations.")
## Polynomial
fn <- function(x) sin(x) + cos(x)
expect_warning(minimaxApprox(fn, -1.5, 1.5, 15L), wrnMess)
## Rational
# The various CRAN and Github testbeds are diverse enough that I cannot find a
# rational minimax approximation example "close enough" to machine precision to
# pass on all of them.
# (AA: 2023-09-01)

# Test consecutive unchanging check and message
fn <- function(x) exp(x) - 1
i <- 5L
opts <- list(conviter = i)
wrnMess <- paste(i, "successive calculated errors were too close to each other",
                 "to warrant further iterations.\n")
## Polynomial
expect_warning(minimaxApprox(fn, -1, 1, 12L, opts = opts), wrnMess)
## Rational
expect_warning(minimaxApprox(fn, -1, 1, c(3L, 3L), opts = opts), wrnMess)

# Test function choosing basis x as 0 trap
wrnMess <- "functional value is 0"
# Polynomial
## Zero is lower bound
expect_warning(minimaxApprox(atan, 0, 1, 14, TRUE), wrnMess)
## Zero is upper bound
fn <- function(x) exp(cos(x)) - 1
expect_warning(minimaxApprox(fn, 0, pi / 2, 4, TRUE), wrnMess)
# Rational
expect_warning(minimaxApprox(sin, 0, pi / 4, c(1L, 1L), TRUE), wrnMess)
## Zero is in the middle. Cheat by using rational where we can pass a known 0.
## For some reason, Github's Mac dies with an error and the Ubuntu/Windows
## servers do not get the Zero basis error. Probably BLAS related, so I will
## only run this at home.
xi <- c(-pi, -2.85, -2.07, -pi / 2, -0.77, -0.2, 0)
if (Sys.info()["nodename"] == "HOME") {
  expect_warning(minimaxApprox(fn, -pi, 0, c(1, 4), TRUE, xi = xi), wrnMess)
}

# Test passing incorrect degree (at minimaxApprox level)
errMsg <- paste("Polynomial approximation takes one value for degree and",
                 "rational approximation takes a vector of two values for",
                 "numerator and denominator degrees. Any other inputs are",
                 "invalid.")
expect_error(minimaxApprox(exp, -1, 1, 1:3), errMsg)

# Test passing xi
## Polynomial - Check that it is ignored
wrnMess <- paste("Polynomial approximation uses Chebyshev nodes for initial",
                 "guess. Any passed xi is ignored.")
expect_message(minimaxApprox(exp, -1, 1, 10L, xi = 6), wrnMess)
## Rational - Check that proper length is passed
errMsg <- paste("Given the requested degrees for numerator and denominator,",
                 "the x-vector needs to have 8 elements.")
xi <- minimaxApprox:::chebNodes(5L, -1, 1)
expect_error(minimaxApprox(exp, -1, 1, c(3L, 3L), xi = xi), errMsg)
# Test that passing proper size works for rational
xi <- xi + 0.01
expect_silent(minimaxApprox(exp, -1, 1, c(2L, 1L), xi = xi))

# Test checkDenom error message
expect_error(minimaxApprox(sin,  0.75 * pi, 1.25 * pi, c(2L, 3L)),
             "The 3 degree polynomial in the denominator has a zero at 2")

## The tests below pass R mac builder AND the Github mac, but for some reason do
## NOT pass CRAN's own mac x86_64 testbed nor on Professor Ripley's Fedora-based
## OpenBLAS platform, so will only run on Windows for now.

if ("windows" %in% tolower(Sys.info()[["sysname"]])) {
  # Test HW Borchers request of returning n degree if n fails but n + 1 works
  # with uppermost effectively 0 with Runge function between -1 and 1 and degree
  # 10.
  ## Test successful restart
  ## Below also tests the failover to QR
  mess <- paste("The algorithm failed while looking for a polynomial of degree",
                "10 but successfully completed when looking for a polynomial",
                "of degree 11 with the largest coefficient's contribution to",
                "the approximation <= the tailtol option. The result is a",
                "polynomial of length 10 as the uppermost coefficient is",
                "effectively 0.")
  fn <- function(x) 1 / (1 + (5 * x) ^ 2)
  control <- c(0.934077073, 0.0, -11.553015692, 0.0, 59.171892231,
               0.0, -134.155250367, 0.0, 135.795965068, 0.0, -50.221129702)
  controlE <- 0.06592293
  expect_message(minimaxApprox(fn, -1, 1, 10L), mess)
  PP <- suppressMessages(minimaxApprox(fn, -1, 1, 10L))
  expect_equal(PP$a, control, tolerance = tol)
  expect_equal(PP$ExpErr, controlE, tolerance = 1e-7) # Only 8 digits in email
  expect_equal(PP$ObsErr, controlE, tolerance = 1e-7) # Only 8 digits in email
}

## Test unsuccessful restart due to two failures
errMsg <- paste("The algorithm neither converged when looking for a",
                 "polynomial of length 22 nor when looking for a polynomial of",
                 "degree 23.")

## Below case has failover to QR
expect_error(minimaxApprox(sin, 0, pi / 2, 22L), errMsg)

# Test tailtol NULL
errMsg <- paste("The algorithm did not converge when looking for a",
                 "polynomial of degree 22 and NULL was passed to the tailtol",
                 "option.")
expect_error(minimaxApprox(sin, 0, pi / 2, 22L, opts = list(tailtol = NULL)),
             errMsg)

## Test unsuccessful restart due to one failures and n + 1 not 0. This must be
## sensitive to precision as it fails on some of github's test platforms, so
## only test on my machine and sacrifice the 100% coverage.
## Below case has failover to QR
if (Sys.info()["nodename"] == "HOME") {
  errMsg <- paste("The algorithm did not converge when looking for a",
                   "polynomial of length 22 and when looking for a polynomial",
                   "of degree 23 the uppermost coefficient is not effectively",
                   "zero.")
  expect_error(minimaxApprox(abs, -0.15, 0.15, 22L), errMsg)
}

# Test ztol
## Polynomial
PP1 <- minimaxApprox(sin, -1, 1, 4L)
PP2 <- minimaxApprox(sin, -1, 1, 4L, opts = list(ztol = 1e-12))

expect_equal(PP2$a[c(2L, 4L)], PP1$a[c(2L, 4L)], tolerance = tol)
expect_identical(PP2$a[c(1L, 3L)], c(0, 0))
expect_equal(PP2$ExpErr, PP1$ExpErr, tolerance = tol)
expect_equal(PP2$ObsErr, PP1$ObsErr, tolerance = tol)
expect_equal(PP2$Basis, PP1$Basis, tolerance = tol)

# This should test RATIONAL failover to QR
expect_error(minimaxApprox(sin, 0, pi / 2, c(100L, 0L)))

################################################################################
# Test minimaxEval
x <- seq(0.1, 0.4, 0.025)
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$ExpErr))
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L))
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$ExpErr))

## Check error trap for mmA object
errMsg <- "This function only works with 'minimaxApprox' objects."
expect_error(minimaxEval(x, sin), errMsg)

## Check not selecting proper basis
errMsg <- "Select either the 'M'onomial or 'C'hebyshev basis."
expect_error(minimaxEval(x, mmA, basis = "A"), errMsg)
expect_error(minimaxEval(x, mmA, basis = 4), errMsg)

## Check asking for Chebyshev when only monomial was run
errMsg <- "Analysis was run using a monomial basis."
expect_error(minimaxEval(x, mmA, basis = "Cheb"), errMsg)

## Check asking for Chebyshev when Chebyshev was run
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L), basis = "c")
expect_true(all(exp(x) - minimaxEval(x, mmA, basis = "c") <= mmA$ExpErr))

## Check asking for monomial when Chebyshev was run
expect_true(all(exp(x) - minimaxEval(x, mmA, basis = "m") <= mmA$ExpErr))

################################################################################
# Test minimaxErr
x <- seq(0.1, 0.4, 0.025)
## Absolute
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_identical(minimaxEval(x, mmA) - exp(x), minimaxErr(x, mmA))

## Relative
mmA <- minimaxApprox(exp, 0, 0.5, 5L, TRUE, basis = "c")
expect_equal((minimaxEval(x, mmA) - exp(x)) / exp(x), minimaxErr(x, mmA),
             tolerance = tol)

## Check error trap
errMsg <- "This function only works with 'minimaxApprox' objects."
expect_error(minimaxErr(x, sin), errMsg)
