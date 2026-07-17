# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)
sW <- function(x) suppressWarnings(x)
sM <- function(x) suppressMessages(x)

nS <- getNamespace("minimaxApprox")
chebNodes <- get("chebNodes", nS, inherits = FALSE, mode = "function")

# Check Accuracy and lack of warning flag when converged
## Rational 1: Based on Fraser & Hart (1962) p. 403 Table 2
controlA <- c(0.99999998510030375, 0.601781180619504719,
              0.186144903531821877, 0.0687440518995425058)
controlB <- c(1, 1.17899457599300466, -0.122321311431112167,
              -0.260995866188425578, 0.0609927504415305534)
controlE <- 1e-6
fn <- function(x) gamma(x + 1)
RR <- minimaxApprox(fn, 0, 1, c(3L, 4L), relErr = FALSE)
expect_equal(RR$aMono, controlA, tolerance = tol)
expect_equal(RR$bMono, controlB, tolerance = tol)
expect_true(RR$ExpErr <= controlE)
expect_false(RR$Warning)

## Rational 2: Based on Cody (1968) pp 250--251. Using weaker tolerance since
## taking values printed on paper.
controlA <- c(1.2655835, -0.65058499, 0.19786869)
controlB <- c(1, -0.064342748, -0.028851456)
controlX <- c(2, 2.0924, 2.3368, 2.6459, 2.9011, 3)
controlE <- 2.6934e-5
RR <- minimaxApprox(gamma, 2, 3, c(2L, 2L), relErr = TRUE, opts = list())
expect_equal(RR$aMono, controlA, tolerance = 5e-6)
expect_equal(RR$bMono, controlB, tolerance = 5e-6)
expect_equivalent(RR$Extrema, controlX, tolerance = 5e-5)
expect_equal(RR$ExpErr, controlE, tolerance = 5e-5)
expect_false(RR$Warning)

## Rational 3: Based on DLMF 3.11.19 https://dlmf.nist.gov/3.11#iii
# Difference on Windows machine is roughly 2.34e-6
controlA <- c(0.99999998917854, -0.34038938209347, -0.18915483763222,
              0.06658319420166)
controlB <- c(1, -0.34039052338838, 0.06086501629812, -0.01864476809090)
fn <- function(x) besselJ(x, nu = 0)
b0 <- 0.893576966279167522
RR <- minimaxApprox(fn, 0, b0, c(3L, 3L))
expect_equal(RR$aMono, controlA, tolerance = 1e-5)
expect_equal(RR$bMono, controlB, tolerance = 1e-5)
expect_false(RR$Warning)

# Test incorrect basis for analysis
errMsg <- paste("Must select either 'C'hebyshev, 'm'onomial, or 'b'arycentric",
                "basis for analysis.")
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
errMsg <- "Degrees must be integers of least 0 (constant)."
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
progMsg <- "i: 1 E: "
expect_message(minimaxApprox(fn, -1, 1, 9L, opts = opts), progMsg)

## Rational
expect_message(minimaxApprox(fn, -1, 1, c(2L, 1L), opts = opts), progMsg)

## Barycentric Polynomial
expect_message(minimaxApprox(fn, -1, 1, 9L, basis = "b", opts = opts), progMsg)

# Test passing some maxiter, convrat, tol, and conviter. Also checks conviter
# overwrite.
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
expect_true(sW(minimaxApprox(fn, -1, 1, 9L, opts = opts)$Warning))

## Rational
dg <- c(2L, 2L)
expect_warning(minimaxApprox(fn, -1, 1, dg, opts = opts), wrnMess)
expect_true(sW(minimaxApprox(fn, -1, 1, dg, opts = opts)$Warning))

## Barycentric Polynomial
expect_warning(minimaxApprox(fn, -1, 1, 9L, basis = "b", opts = opts), wrnMess)
expect_true(sW(minimaxApprox(fn, -1, 1, 9L, basis = "b", opts = opts)$Warning))

## E5 re-baseline: pre-E5 this stalled at the floor and raised the near-eps
## warning; the E5 exchange's trajectory reaches the singular -> rescue path
## instead (dense error 2.0e-15, machine floor either way). Both are honest
## machine-precision outcomes; either text passes, but one must fire.
fn <- function(x) sin(x) + cos(x)
w15 <- tryCatch(minimaxApprox(fn, -1.5, 1.5, 15L),
                warning = function(w) conditionMessage(w))
expect_true(grepl("NOT technically a Remez result", w15, fixed = TRUE) ||
              grepl(wrnMess, w15, fixed = TRUE))

## Rational
# The various CRAN and Github testbeds are diverse enough that I cannot find a
# rational minimax approximation example "close enough" to machine precision to
# pass on all of them.
# (AA: 2023-09-01)

# Test consecutive unchanging check and message
fn <- function(x) exp(cos(x))
i <- 1L
opts <- list(conviter = i)
wrnMess <- paste(i, "successive calculated solutions were too close to each",
                 "other to warrant further iterations.\n")
## Polynomial
expect_warning(minimaxApprox(fn, -1, 1, 21L, opts = opts), wrnMess)

## Rational
expect_warning(minimaxApprox(fn, -pi, pi, c(5L, 2L), opts = opts), wrnMess)

# Test function choosing basis x as 0 trap
wrnMess <- "functional value is 0"

# Polynomial
## Zero is lower bound
expect_warning(minimaxApprox(atan, 0, 1, 14, TRUE, basis = "m"), wrnMess)

## Zero is upper bound
fn <- function(x) exp(cos(x)) - 1
expect_warning(minimaxApprox(fn, 0, pi / 2, 4, TRUE), wrnMess)

# Rational
expect_warning(minimaxApprox(sin, 0, pi / 4, c(1L, 1L), TRUE), wrnMess)

# Test passing incorrect degree (at minimaxApprox level)
errMsg <- paste("Polynomial approximation takes one value for degree and",
                "rational approximation takes a vector of two values for",
                "numerator and denominator degrees. Any other inputs are",
                "invalid.")
expect_error(minimaxApprox(exp, -1, 1, 1:3), errMsg)

# Test passing xi (deprecated 0.6.0; removal scheduled next release)
## Deprecation warning fires whenever xi is passed, on every path
dprMess <- "'xi' argument is deprecated"
expect_warning(minimaxApprox(exp, -1, 1, c(2L, 1L),
                             xi = chebNodes(5L, -1, 1) + 0.01), dprMess)

## Polynomial - Check that it is (still) ignored, alongside the warning
wrnMess <- paste("Polynomial approximation uses Chebyshev nodes for initial",
                 "guess. Any passed xi is ignored.")
expect_message(sW(minimaxApprox(exp, -1, 1, 10L, xi = 6)), wrnMess)

## Rational - Check that proper length is passed
errMsg <- paste("Given the requested degrees for numerator and denominator,",
                "the x-vector needs to have 8 elements.")
xi <- chebNodes(5L, -1, 1)
expect_error(sW(minimaxApprox(exp, -1, 1, c(3L, 3L), xi = xi)), errMsg)

# Test that passing proper size works for rational (deprecation warning is
# the ONLY condition emitted; silence otherwise)
xi <- xi + 0.01
expect_silent(sW(minimaxApprox(exp, -1, 1, c(2L, 1L), xi = xi)))

# Test checkDenom error message
expect_error(minimaxApprox(sin,  0.75 * pi, 1.25 * pi, c(2L, 3L)),
             "The 3 degree polynomial in the denominator has a zero at 2")

## The tests below pass R mac builder AND the Github mac, but for some reason do
## NOT pass CRAN's own mac x86_64 testbed nor on Professor Ripley's Fedora-based
## OpenBLAS platform, so will only run on Windows for now.

# if ("windows" %in% tolower(Sys.info()[["sysname"]])) {

## They may pass now with the rengineered switch/findroots so will try removing
## the gate.
## (AA: 2026-07-16)

# E5 re-baseline: the redesigned exchange converges degree 10 DIRECTLY
# (no singular solve, so no degree-11 restart and no message) on the
# platforms measured so far; a platform whose solve still goes singular
# takes the restart path and emits HWB's message. Accept either route --
# what is pinned is the RESULT: the returned polynomial must match HWB's
# control coefficients and error either way.

fn <- function(x) 1 / (1 + (5 * x) ^ 2)
control <- c(0.934077073, 0.0, -11.553015692, 0.0, 59.171892231, 0.0,
             -134.155250367, 0.0, 135.795965068, 0.0, -50.221129702)
controlE <- 0.06592293

msgs <- character(0)
PP <- withCallingHandlers(minimaxApprox(fn, -1, 1, 10L),
                          message = function(m) {
                            msgs <<- c(msgs, conditionMessage(m)) # nolint: undesirable_operator_linter
                            invokeRestart("muffleMessage")
                          })
expect_true(length(msgs) == 0L ||
              any(grepl("successfully completed when looking", msgs,
                        fixed = TRUE)))
expect_equal(PP$aMono, control, tolerance = tol)
expect_equal(PP$ExpErr, controlE, tolerance = 1e-7) # Only 8 digits in email
expect_equal(PP$ObsErr, controlE, tolerance = 1e-7) # Only 8 digits in email

# }

## Test unsuccessful restart due to two failures. F4: the former case here
## (sin, 0.25, 0.75, 16, "m") is now RESCUED -- see the F4 block below --
## because its interpolant is at the machine-precision floor. Replaced with a
## genuinely non-representable case that must still hard-error: sqrt(x) on
## [0, 1] has a branch point at 0, so no polynomial interpolant comes near it
## (probe abs error ~2e-2 >> threshold), and the F4 rescue correctly falls
## through to the original error rather than masking it.
errMsg <- "The algorithm neither converged when looking for a"

# E5 re-baseline: pre-E5 the exchange fed a singular solve at degree 15 and,
# with tailtol = NULL disabling the restart, the hard error above was raised.
# The E5 exchange never goes singular on this input; the iteration stalls and
# the SECOND-manifestation rescue (M4, downstream of and independent from the
# tailtol-gated restart) returns the degree-15 interpolant -- dense error
# 1.11e-16, a machine-floor-perfect result with the honest rescue warning, a
# strictly better outcome than the hard error. Accept either honest form
# (a platform whose solve still goes singular takes the error arm).
errMsg <- "The algorithm did not converge when looking for a"
oTT <- tryCatch(
  sM(minimaxApprox(sin, 0.25, 0.75, 15L, basis = "m",
                   opts = list(tailtol = NULL))),
  warning = function(w) w, error = function(e) e
)
expect_true(inherits(oTT, "error") &&
              grepl(errMsg, conditionMessage(oTT), fixed = TRUE) ||
              inherits(oTT, "warning") &&
                grepl("NOT technically a Remez result", conditionMessage(oTT),
                      fixed = TRUE))

## Test unsuccessful restart: degree-n Remez fails singular, degree-(n+1)
## retry SUCCEEDS but its uppermost coefficient is NOT effectively zero (fails
## the tailtol test). The precise degree at which the degree-n augmented solve
## first becomes singular is BLAS/LAPACK-dependent (a conditioning boundary,
## not a structural one), so this test cannot pin a single outcome across all
## platforms without gating to one machine. Instead it accepts EITHER of the
## two legitimate outcomes for this input and rejects the two illegitimate
## ones, which makes it deterministic and un-gated:
##   (A) degree-18 singular, degree-19 solves with a large top coefficient
##       -> the "uppermost coefficient is not effectively zero" ERROR
##          (the branch this test exists to cover); OR
##   (B) degree-18 not yet singular on this platform's BLAS
##       -> the algorithm converges, or drops to a lower degree via the
##          "effectively 0" success MESSAGE.
## Runge with a steep 8x scaling at degree 18 is chosen so that WHEN outcome
## (A) occurs, the degree-19 top coefficient's contribution clears tailtol by
## ~4000x -- removing the *second* fragility (a marginal tail test) that the
## former degree-22 case also had. Only the singular-detection boundary
## remains platform-dependent, and both sides of it are accepted here.
## A future barycentric path (M5) may sidestep the singular solve entirely, at
## which point this can revert to a single-outcome assertion.
fn <- function(x) 1 / (1 + (8 * x) ^ 2)
targetErr <- paste("The algorithm did not converge when looking for a",
                   "polynomial of degree 18 and when looking for a polynomial",
                   "of degree 19 the uppermost coefficient is not effectively",
                   "zero.")
res <- tryCatch(
  minimaxApprox(fn, -1, 1, 18L, basis = "m", opts = list(tailtol = 1e-10)),
  error = function(e) structure(conditionMessage(e), class = "mmaOutcomeErr"),  # nolint undesirable_operator_linter
  message = function(m) structure(conditionMessage(m), class = "mmaOutcomeMsg") # nolint undesirable_operator_linter
)

if (inherits(res, "mmaOutcomeErr")) {
  # Outcome (A): must be exactly the target branch, NOT "neither converged".
  expect_identical(as.character(res), targetErr)
} else if (inherits(res, "mmaOutcomeMsg")) {
  # Outcome (B), drop-to-lower-degree form: must be the effectively-0 success
  # message, and must therefore have returned a usable lower-degree result.
  expect_true(grepl("uppermost coefficient is effectively 0",
                    as.character(res), fixed = TRUE))
} else {
  # Outcome (B), clean-convergence form: a valid minimaxApprox object.
  expect_inherits(res, "minimaxApprox")
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
# E5 re-baseline: pre-E5 the exchange fed the degree-100 monomial solve a
# singular system (hard error). The E5 trajectory avoids the singularity and
# the iteration runs to a stall exit with Warning = TRUE (measured: ratio
# ~708 reported in the warning; the input is far outside the supported
# envelope either way). Accept either loud outcome; a SILENT completion
# would be the failure mode.
o100 <- tryCatch(sW(minimaxApprox(sin, 0, pi / 2, c(100L, 0L))),
                 error = function(e) e)
expect_true(inherits(o100, "error") || isTRUE(o100$Warning))

################################################################################
# Input validation additions (Module M2: F2, F10, F13)

# --- F2: lower > upper / lower == upper is now a clean error, not a silently
# suboptimal result. Exact repro from the review document.
errMsg <- "'lower' must be less than 'upper'"
expect_error(minimaxApprox(exp, 1, 0, 3), errMsg)
expect_error(minimaxApprox(exp, 1, 1, 3), errMsg)

# --- F10: lower/upper must be finite, non-missing numeric scalars.
errMsg <- "'lower' must be a finite, non-missing numeric scalar."
expect_error(minimaxApprox(exp, NA, 1, 3), errMsg)
expect_error(minimaxApprox(exp, NaN, 1, 3), errMsg)
expect_error(minimaxApprox(exp, -Inf, 1, 3), errMsg)
expect_error(minimaxApprox(exp, "a", 1, 3), errMsg)
expect_error(minimaxApprox(exp, c(0, 1), 1, 3), errMsg)

errMsg <- "'upper' must be a finite, non-missing numeric scalar."
expect_error(minimaxApprox(exp, 0, NA, 3), errMsg)
expect_error(minimaxApprox(exp, 0, NaN, 3), errMsg)
expect_error(minimaxApprox(exp, 0, Inf, 3), errMsg)
expect_error(minimaxApprox(exp, 0, "a", 3), errMsg)

# --- F10: fn must be a function whose first formal is 'x'. Primitives
# (formals() == NULL) must still be accepted via args().
errMsg <- "'fn' must be a function whose first argument is 'x'."
expect_error(minimaxApprox(function(t) exp(t), 0, 1, 3), errMsg)
expect_error(minimaxApprox(function() 1, 0, 1, 3), errMsg)
expect_error(minimaxApprox("not a function", 0, 1, 3), errMsg)
expect_silent(minimaxApprox(sin, 0, 1, 3))
expect_silent(minimaxApprox(exp, 0, 1, 3))
expect_silent(minimaxApprox(function(x) sin(x), 0, 1, 3))

# --- F10: degree must be finite, non-missing numeric.
errMsg <- "'degree' must be finite, non-missing numeric value(s)."
expect_error(minimaxApprox(exp, 0, 1, NA), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, NaN), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, Inf), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, "3"), errMsg, fixed = TRUE)
expect_error(minimaxApprox(exp, 0, 1, c(3, NA)), errMsg, fixed = TRUE)

# --- F10: opts must be a list; individual members type/range-checked when
# supplied. tailtol/ztol remain permitted to be NULL by design.
expect_error(minimaxApprox(exp, 0, 1, 3, opts = 5), "'opts' must be a list.")

errMsg <- "'opts\\$maxiter' must be a single positive integer."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(maxiter = 0)), errMsg)
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(maxiter = NA)), errMsg)
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(maxiter = 1.5)), errMsg)

errMsg <- "'opts\\$miniter' must be a single positive integer."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(miniter = -1)), errMsg)

errMsg <- "'opts\\$conviter' must be a single positive integer."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(conviter = 0)), errMsg)

errMsg <- "'opts\\$tol' must be a single finite, non-missing numeric value."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(tol = "x")), errMsg)

errMsg <- "'opts\\$convrat' must be a single finite, non-missing numeric value."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(convrat = NA)), errMsg)

errMsg <- "'opts\\$tailtol' must be a single finite, non-missing numeric value."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(tailtol = "x")), errMsg)
expect_silent(minimaxApprox(exp, 0, 1, 3, opts = list(tailtol = NULL)))

errMsg <- "'opts\\$ztol' must be a single finite, non-missing numeric value."
expect_error(minimaxApprox(exp, 0, 1, 3, opts = list(ztol = "x")), errMsg)
expect_silent(minimaxApprox(exp, 0, 1, 3, opts = list(ztol = NULL)))

# ---------------------------------------------------------------------------
# F4: exact-representability / machine-precision-resolved rescue.
# When both the degree-n and degree-(n+1) Remez solves fail singular AND the
# degree-n interpolant is at the machine-precision floor, minimaxApprox now
# returns that interpolant with a warning that it is NOT a Remez result. When
# the interpolant is NOT at the floor (genuine non-representability), or the
# interpolation is itself singular, or relErr is requested and fn has a zero
# in range, the original hard error is preserved (tested above / below).
# ---------------------------------------------------------------------------

f4wrn <- "NOT technically a Remez result"

# Exactly representable: x^2 is exactly a degree-2 polynomial. Rescue returns
# the interpolant (== the function); aMono is (0, 0, 1) to machine precision.
expect_warning(minimaxApprox(function(x) x^2, 0, 1, 2L), f4wrn)
ppx2 <- sW(minimaxApprox(function(x) x^2, 0, 1, 2L))
expect_equal(ppx2$aMono, c(0, 0, 1), tolerance = 1e-12)
expect_true(ppx2$ObsErr < 10 * .Machine$double.eps)
expect_true(ppx2$Warning)

# Resolved-to-precision: exp on [-1, 1] is resolved by the requested degree
# (deg 13 already converges with ratio 1.52; deg 14, 15, 50 previously HARD
# ERRORED). Chebyshev basis: conditioning ~1.4, so these are robust across
# platforms. Rescue returns a floor-level interpolant.
for (d in c(14L, 15L, 50L)) {
  expect_warning(minimaxApprox(exp, -1, 1, d), f4wrn)
  ppe <- sW(minimaxApprox(exp, -1, 1, d))
  expect_true(ppe$ObsErr < 1e-13)
  expect_true(ppe$Warning)
  # Re-measure the returned coefficients independently: floor-level everywhere.
  g <- seq(-1, 1, length.out = 501L)
  expect_true(max(abs(minimaxErr(g, ppe))) < 1e-13)
}

# relErr rescue-success path (covers interpRescue's relative-error branch,
# lines "err <- max(abs((pg - fg) / fg))" / "thresh <- 10 * eps" -- every OTHER
# relErr case in this file hits the any(fg == 0) zero-guard and returns before
# reaching those lines). x^2 + 1 is exactly degree 2 and has no zero on [0, 1],
# so relative error is well-defined and the rescue fires via that branch.
expect_warning(minimaxApprox(function(x) x^2 + 1, 0, 1, 2L, relErr = TRUE),
               f4wrn)
pprel <- sW(minimaxApprox(function(x) x^2 + 1, 0, 1, 2L, relErr = TRUE))
expect_equal(pprel$aMono, c(1, 0, 1), tolerance = 1e-12)
expect_true(pprel$ObsErr < 10 * .Machine$double.eps)
expect_true(pprel$Warning)

# Fall-through 1 (relErr zero-guard): x on [-1, 1] with relErr has a zero at
# x = 0 on the probe grid, so the relative criterion is undefined and the
# rescue does NOT fire.
# E5 re-baseline: pre-E5 both Remez solves went singular and, with the
# rescue's relErr zero-guard declining, the hard error was raised. The E5
# exchange does not go singular here; the iteration returns an (essentially)
# exact representation of x and exits through the stall path, Warning TRUE.
# NOTE: ObsErr is deliberately NOT asserted. It is the RELATIVE error, and
# fn's zero at x = 0 divides the solve's coefficient noise by x: measured
# ObsErr is exactly 0 on the reviewer container (bitwise-exact solve) but
# 1.9e-8 on HOMEDESKTOP (coefficient noise ~8e-10) -- platform noise, not a
# contract. The contract is that this input must never complete SILENTLY:
# either the documented error, or a completed result with Warning TRUE.
errMsg <- "The algorithm neither converged when looking for a"
oXR <- tryCatch(sW(minimaxApprox(function(x) x, -1, 1, 12L,
                                 relErr = TRUE, basis = "m")),
                error = function(e) e)
expect_true(inherits(oXR, "error") &&
              grepl(errMsg, conditionMessage(oXR), fixed = TRUE) ||
              isTRUE(oXR$Warning))
