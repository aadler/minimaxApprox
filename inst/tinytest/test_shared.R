# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)
sW <- function(x) suppressWarnings(x)
sM <- function(x) suppressMessages(x)

nS <- getNamespace("minimaxApprox")
fC <- get("fC", nS, inherits = FALSE, mode = "function")
chebNodes <- get("chebNodes", nS, inherits = FALSE, mode = "function")
callFun <- get("callFun", nS, inherits = FALSE, mode = "function")
isOscil <- get("isOscil", nS, inherits = FALSE, mode = "function")
isConverged <- get("isConverged", nS, inherits = FALSE, mode = "function")
zb <- get("zeroBasisPerturb", nS, inherits = FALSE, mode = "function")
evalFunc <- get("evalFunc", nS, inherits = FALSE, mode = "function")
remPoly <- get("remPoly", nS, inherits = FALSE, mode = "function")
remRat <- get("remRat", nS, inherits = FALSE, mode = "function")
remErr <- get("remErr", nS, inherits = FALSE, mode = "function")
polyCoeffs <- get("polyCoeffs", nS, inherits = FALSE, mode = "function")
findRoots <- get("findRoots", nS, inherits = FALSE, mode = "function")
ratCoeffs <- get("ratCoeffs", nS, inherits = FALSE, mode = "function")
switchX <- get("switchX", nS, inherits = FALSE, mode = "function")
checkDenom <- get("checkDenom", nS, inherits = FALSE, mode = "function")
isUnchanging <- get("isUnchanging", nS, inherits = FALSE, mode = "function")
checkIrrelevant <- get("checkIrrelevant", nS, inherits = FALSE,
                       mode = "function")
tailContribution <- get("tailContribution", nS, inherits = FALSE,
                        mode = "function")
rlc <- get("refLocalCheck", nS, inherits = FALSE, mode = "function")
REFLOCALTOL <- get("REFLOCALTOL", nS, inherits = FALSE, mode = "double")

opts <- list(maxiter = 100L, miniter = 10L, conviter = 10L,
             showProgress = FALSE, convrat = 1.000000001, tol = 1e-14,
             ztol = .Machine$double.eps)

# Test fC
expect_identical(fC(1.234567, f = "e"), "1.234567e+00")
expect_identical(fC(1.234567, d = 2, f = "e"), "1.23e+00")

# Test chebNodes
n <- 6L
k <- seq_len(n) - 1L
# See https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema
control <- sort(cos(pi * (k + 0.5) / n))
expect_equal(chebNodes(n, -1, 1), control, tolerance = tol)
expect_equal(chebNodes(6.2, -1, 1), chebNodes(n, -1, 1), tolerance = tol)

# Test callFun
## Test functionality
fn <- function(x) tan(x) - x ^ 3
control <- tan(-0.4) - (-0.4) ^ 3
expect_equal(callFun(fn, -0.4), control, tolerance = tol)

## Test error trapping
expect_error(callFun("x ^ 2", -0.4), "Unable to parse function.")

# Test isOscil
control <- c(-2, 1, -3, 4, -1, 6, -7)
expect_true(isOscil(control))

control <- c(-2, 1, -3, 4, -1, -6)
expect_false(isOscil(control))

# Test evalFunc
# Tests using Chebyshev polynomials are currently in test_Chebyshev.R
x <- c(-0.1, 0.2, 2)
controlN <- 1 + 2 * x + 3 * x ^ 2 + 4 * x ^ 3

## Polynomial
P <- list(a = 1:4)
expect_equal(evalFunc(x, P, "m", -1, 1), controlN, tolerance = tol)

## Rational
R <- list(a = 1:4, b = c(1, 2.2, 4.1))
controlD <- 1 + 2.2 * x + 4.1 * x ^ 2
control <- controlN / controlD
expect_equal(evalFunc(x, R, "m", -1, 1), control, tolerance = tol)

# Test remErr
# Using fact that exp(1) has analytic answer for degree 1 and pass a zero-degree
# polynomial in the denominator for the rational test
fn <- function(x) exp(x)
m <- exp(1) - 1
c <- (exp(1) - m * log(m)) / 2
tstFn <- function(x) m * x + c
x <- chebNodes(3, 0, 1)
control <- tstFn(x) - exp(x)

## Polynomial
PP <- remPoly(fn, 0, 1, 1, FALSE, "m", opts)
expect_equal(remErr(x, PP, fn, FALSE, "m", 0, 1), control, tolerance = tol)
## Rational
RR <- remRat(fn, 0, 1, 1, 0, FALSE, "m", NULL, opts)
expect_equal(remErr(x, RR, fn, FALSE, "m", 0, 1), control, tolerance = tol)

# Test findRoots
## This one will rely on expm1(x) and exp(x) - 1 being close
fn <- function(x) exp(x) - 1
x <- chebNodes(3, 0, 1)

## Polynomial
QQ <- polyCoeffs(x, expm1, TRUE, "m", 0, 1, opts$ztol)
control <- findRoots(x, QQ, expm1, TRUE, "m", 0, 1)
PP <- polyCoeffs(x, fn, TRUE, "m", 0, 1, opts$ztol)
r <- findRoots(x, PP, fn, TRUE, "m", 0, 1)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1e-7)

## Rational
QQ <- ratCoeffs(x, 0, expm1, 1L, 0L, TRUE, "m", 0, 1, opts$ztol)
control <- findRoots(x, QQ, expm1, TRUE, "m", 0, 1)
RR <- ratCoeffs(x, 0, fn, 1L, 0L, TRUE, "m", 0, 1, opts$ztol)
r <- findRoots(x, RR, fn, TRUE, "m", 0, 1)
## Need weaker tolerance here since functions are not exactly the same
expect_equal(r, control, tolerance = 1e-7)

## E5 contract tests. (The pre-E5 "error trap" test here passed an UNDEFINED
## object A: uniroot errored with "object not found", the old tryCatch
## swallowed it, and the fallback path returned 1.2 -- a latent test bug that
## only tested the trap by accident. The fallback no longer exists.)
## A rootless error curve now yields numeric(0), not a substituted endpoint.
pc <- list(a = 5)
fc <- function(x) 3
expect_identical(findRoots(c(0.4, 0.6), pc, fc, FALSE, "m", 0, 1), double(0))
## D1 regression (MP2 B1-1): a sign change OUTSIDE the between-reference
## span -- here the root at 0.2 lies left of the reference {0.5, 0.6} -- was
## structurally invisible pre-E5 and must now be found.
pl <- list(a = c(-0.2, 1))
f0 <- function(x) 0 * x
r <- findRoots(c(0.5, 0.6), pl, f0, FALSE, "m", 0, 1)
expect_equal(r, 0.2, tolerance = 1e-7)

# Test switchX
# Assuming function is correct, replicate a previous result.
## Polynomial
control <- c(-1, 0.10264791208519766, 0.33735881337846646, 0.62760501759598053,
             0.88066205512236839, 1)
fn <- function(x) sin(x) + cos(x)
x <- chebNodes(6, 0, 1)
PP <- polyCoeffs(x, fn, FALSE, "m", 0, 1, opts$ztol)
r <- findRoots(x, PP, fn, FALSE, "m", 0, 1)
x <- switchX(r, -1, 1, PP, fn, FALSE, "m", x)
# Need weaker tolerance here due to different build platforms
expect_equivalent(x, control, tolerance = 3.5e-5)

## Rational
control <- c(-1, -0.6706726462230721, -2.8931353340360859e-14,
             0.67067262060160282, 1)
fn <- function(x) ifelse(abs(x) < 1e-20, 1, sin(x) / x)
x <- chebNodes(5, -1, 1)
RR <- ratCoeffs(x, 0, fn, 2L, 1L, FALSE, "m", -1, 1, opts$ztol)
r <- findRoots(x, RR, fn, FALSE, "m", -1, 1)
x <- switchX(r, -1, 1, RR, fn, FALSE, "m", x)
# Need weaker tolerance here due to different build platforms
expect_equivalent(x, control, tolerance = 3.5e-5)

## E5: degenerate (too-few-candidates) fallback. A constant error curve
## offers a single sign region -- fewer alternating candidates than the
## reference size -- so switchX must return the PREVIOUS reference unchanged
## (routing the loop to the existing isUnchanging stall exit; deliberately
## no new degenerate handling, M3 lesson) with the ZeroBasis attribute set.
R <- list(a = 0, b = 1)
fn <- function(x) 3
xk <- c(0.25, 0.75)
xf <- switchX(numeric(0), 0, 1, R, fn, FALSE, "m", xk)
expect_equivalent(as.vector(xf), xk, tolerance = tol)
expect_false(attr(xf, "ZeroBasis"))

fn <- function(x) -3
xf <- switchX(numeric(0), 0, 1, R, fn, FALSE, "m", xk)
expect_equivalent(as.vector(xf), xk, tolerance = tol)
expect_false(attr(xf, "ZeroBasis"))

# Check isConverged
errs <- c(-0.1, 0.1, -0.1)
E <- 0.1
expect_true(isConverged(errs, E, 1.05, 1e-12))
E <- 0.05
expect_false(isConverged(errs, E, 1.05, 1e-12))
E <- 0.1
errs <- c(-0.2, 0.1, -0.1)
expect_false(isConverged(errs, E, 1.05, 1e-12))

# Test checkDenom
# NOTE (M6): these previously passed basis = TRUE, which only ever "worked"
# because switch(EXPR = TRUE, m = ..., ...) coerces TRUE to integer 1 and
# picks the first listed alternative POSITIONALLY, regardless of its name --
# an accident of switch()'s non-character-EXPR behavior, not a valid basis
# value. checkDenom's now-explicit if (basis == "m") dispatch (needed to
# route Chebyshev through the M6 affine map) correctly stops honoring that
# accident. The polynomial -0.5 + x (root 0.5 on [0,1], no root on [1,2]) is
# unambiguously the intended monomial-basis case; fixed to basis = "m".
expect_equal(checkDenom(c(-0.5, 1), 0, 1, "m"), 0.5)
expect_null(checkDenom(c(-0.5, 1), 1, 2, "m"))

# --------------------------------------------------------------------------
# M3 additions: F7, F8, F9, F6, F3
# --------------------------------------------------------------------------

# F7 -- isUnchanging must not flag rapidly-improving errors as stagnation.
errs_last <- rep(1e-3, 4)
convrat <- 1.000000001
tolF <- 1e-14
## 10x-shrink: genuine improvement, must NOT be flagged unchanging.
expect_false(isUnchanging(errs_last / 10, errs_last, convrat, tolF))
## Static vector: genuinely unchanging, must be flagged.
expect_true(isUnchanging(errs_last, errs_last, convrat, tolF))
## Straddling the two-sided band (one ratio far below 1/convrat): must NOT be
## flagged, since not all elements are close to unchanged.
straddle <- c(1e-3, 1e-3, 1e-3, 5e-4)
expect_false(isUnchanging(straddle, errs_last, convrat, tolF))
## Zero-denominator perturbation path still works (both become 1e-12, ratio
## exactly 1, difference exactly 0).
expect_true(isUnchanging(rep(0, 4), rep(0, 4), convrat, tolF))

# F8 -- isOscil must not propagate NA from NaN/NA input.
expect_false(isOscil(c(1, NaN, -1)))
expect_false(isOscil(c(1, NA, -1)))
## Zero error treated as non-oscillating by design.
expect_false(isOscil(c(1, 0, -1)))
## Genuinely alternating case unaffected.
expect_true(isOscil(c(-2, 1, -3, 4)))
## isConverged no longer errors inside an if() when isOscil returns FALSE
## instead of NA.
expect_false({
  ok <- TRUE
  tryCatch(
    if (isConverged(c(1, NaN, -1), 1, 1.05, 1e-12)) NULL,
    error = function(e) ok <<- FALSE) # nolint: undesirable_operator_linter
  !ok
})

# F9 -- zeroBasisPerturb: absolute step at ordinary |x| (released behavior),
# escalating to a magnitude-scaled step only where the absolute step is
# absorbed (large |x|). Tested at both endpoints, both signs, and interior.

## Ordinary |x_i|: interior perturbation is byte-identical to the released
## absolute 1e-12 nudge (regression guard for the machine-precision
## ZeroBasis end-to-end cases). fn = x^2-4 has its zero at the interior
## point -2; released code probed c(-2-1e-12, -2+1e-12) and picked by fn.
oldInterior <- {
  cand <- c(-2 - 1e-12, -2 + 1e-12)
  cand[which.max(cand ^ 2 - 4)]
}
expect_identical(zb(-2, -3, -1, function(x) x ^ 2 - 4, TRUE), oldInterior)
## Same interior case, minimize direction (maximize = FALSE): exercises the
## which.min(fnreplace) branch, otherwise only reachable via the platform-
## sensitive x^2-4 relErr end-to-end case (see the HOMEDESKTOP-gated block in
## test_MiniMaxApprox.R). Deterministic, no Remez iteration/BLAS involved --
## calls zeroBasisPerturb directly.
oldInteriorMin <- {
  cand <- c(-2 - 1e-12, -2 + 1e-12)
  cand[which.min(cand ^ 2 - 4)]
}
expect_identical(zb(-2, -3, -1, function(x) x ^ 2 - 4, FALSE), oldInteriorMin)

## Large |x|: the plain absolute step is a no-op, so the escalated step must
## actually move the point.
x5 <- 5e4
expect_true((x5 - 1e-12) == x5)
expect_true((x5 + 1e-12) == x5)   # confirms the no-op
expect_true(zb(x5, -1e5, 1e6, function(x) x, TRUE) != x5)
## Lower endpoint, both signs: perturbed point strictly inside (l, u).
expect_true(zb(5e4, 5e4, 1e6, function(x) x, TRUE) > 5e4)
expect_true(zb(-5e4, -5e4, 1e5, function(x) x, TRUE) > -5e4)
## Upper endpoint, both signs: perturbed point strictly inside (l, u).
expect_true(zb(5e4, -1e6, 5e4, function(x) x, TRUE) < 5e4)
expect_true(zb(-5e4, -1e5, -5e4, function(x) x, TRUE) < -5e4)
## Ordinary-magnitude endpoints: byte-identical to the released absolute step.
expect_identical(zb(-3, -3, -1, function(x) x, TRUE), -3 + 1e-12)
expect_identical(zb(-1, -3, -1, function(x) x, TRUE), -1 - 1e-12)
## l == 0 / u == 0: pure absolute step, strictly inside.
expect_identical(zb(0, 0, 1, function(x) x, TRUE), 0 + 1e-12)
expect_identical(zb(0, -1, 0, function(x) x, TRUE), 0 - 1e-12)

# F6 -- checkIrrelevant must be basis-aware.
## Confirmed repro: genuine ~1e-3 Chebyshev contributions on [0, 0.5] must
## survive ztol = 1e-4 (previously zeroed by the monomial xmax^k scaling).
aCheb <- c(1, rep(1e-3, 6))
rCheb <- checkIrrelevant(aCheb, 0, 0.5, 1e-4, "c")
expect_equal(rCheb, aCheb, tolerance = tol)
## A genuinely negligible coefficient must still be zeroed, in both bases.
aNeg <- c(1, rep(1e-10, 6))
expect_equal(checkIrrelevant(aNeg, 0, 0.5, 1e-4, "c"),
             c(1, rep(0, 6)), tolerance = tol)
expect_equal(checkIrrelevant(aNeg, 0, 0.5, 1e-4, "m"),
             c(1, rep(0, 6)), tolerance = tol)
## Monomial-basis results bitwise-unchanged versus the pre-F6 formula.
set.seed(20260710)
aRand <- rnorm(8)
oldMonomial <- {
  nn <- length(aRand)
  xmax <- max(abs(0), abs(2))
  ifelse(abs(aRand * xmax ^ (seq_len(nn) - 1L)) <= 1e-6, 0, aRand)
}
expect_identical(checkIrrelevant(aRand, 0, 2, 1e-6, "m"), oldMonomial)

# F3 -- tailContribution: abs() + basis-correct scale.
## A large NEGATIVE top coefficient must now exceed tailtol, in both bases
## (previously silently passed with no abs()).
expect_true(tailContribution(-1e-5, 12, -1, 1, "m") > 1e-10)
expect_true(tailContribution(-1e-5, 12, -1, 1, "c") > 1e-10)
## A tiny coefficient of either sign must not exceed tailtol.
expect_false(tailContribution(-1e-20, 12, -1, 1, "m") > 1e-10)
expect_false(tailContribution(1e-20, 12, -1, 1, "c") > 1e-10)

# CP-1 (E5 Phase 1) -- reference-local certificate on all paths.
# Mechanism: MP2 B1-1 / E5 brief section 0.1 (D3 certification vacuum).

## Unit tests of refLocalCheck itself.
## Healthy classical fit: converged exp deg 6 -- certificate passes.
hf <- sW(minimaxApprox(exp, -1, 1, 6L))
expect_false(rlc(list(a = hf$a), exp, FALSE, "c", -1, 1, hf$ExpErr)$refLocal)
## Synthetic bad fit: a deliberately-wrong polynomial with a small claimed
## leveled error must fail its certificate.
expect_true(rlc(list(a = c(1, 0, 0)), exp, FALSE, "m", -1, 1, 1e-6)$refLocal)
## relErr zero-on-probe guard: sin has an exact zero at the probe point 0
## (2001 points on [-1, 1] include 0), so the check must skip, not fire.
zg <- rlc(list(a = c(1, 0, 0)), sin, TRUE, "m", -1, 1, 1e-6)
expect_false(zg$refLocal)
expect_true(is.na(zg$gridSup))
## Floor gate: an expe below 100 * eps * max(1, ||f||) skips the check even
## for a wrong polynomial (the ratio would be noise at that magnitude).
fg <- rlc(list(a = c(1, 0, 0)), exp, FALSE, "m", -1, 1, 1e-15)
expect_false(fg$refLocal)
expect_true(is.na(fg$gridSup))

## Absolute-excess condition: exp deg 11 Chebyshev converges with a grid
## ratio marginally over REFLOCALTOL (measured 1.0011) whose absolute excess
## is ~4 * eps * ||f|| -- solve noise, not a basin miss. Must stay silent.
o11 <- minimaxApprox(exp, -1, 1, 11L, basis = "c")
expect_false(o11$Warning)

## The B1-1 quintet, post-E5 (exchange redesign): every case must now
## converge WARNING-FREE to the true minimax with a dense-grid certificate,
## on every platform. Oracles: atan[0, 3] deg 3 = 4.802475e-3 (dual-path,
## MP2); atan deg 9 = 1.143854e-5 and sin deg 9 = 2.396019e-11 (HOMEDESKTOP
## good-basin measurements, pre-E5; the sin value = 2 * J_11(1)); atan
## deg 11 = 1.662360e-6 and cos deg 4 = 4.187752e-5 (parity-staircase
## oracles: the degree-(n+1) classical fits agree to displayed digits,
## E_n = E_{n+1} for odd/even functions). A refLocal warning here is a
## certificate failure -- post-E5 that is a bug alarm, not an accepted
## outcome.
certAt <- function(fn, l, u, d, b, oracle) {
  o <- sM(minimaxApprox(fn, l, u, d, basis = b))
  g <- seq(l, u, length.out = 2e5L)
  gridRatio <- max(abs(sM(minimaxEval(g, o) - fn(g)))) / o$ExpErr
  !o$Warning && gridRatio <= REFLOCALTOL && abs(o$ExpErr / oracle - 1) < 1e-4
}
expect_true(certAt(atan, 0, 3, 3L, "b", 4.802475e-3))
expect_true(certAt(atan, -1, 1, 9L, "m", 1.143854e-5))
expect_true(certAt(atan, -1, 1, 11L, "b", 1.662360e-6))
expect_true(certAt(sin, -1, 1, 9L, "c", 2.396019e-11))
## The fifth (CP-1-discovered) case: even fn at even degree, symmetric
## interval -- exercises both the D1/D2 exchange fixes and the E5
## endpoint-candidate rule that breaks the parity-degenerate (h ~ 0)
## symmetric-reference fixed point.
expect_true(certAt(cos, -1, 1, 4L, "b", 4.187752e-5))

### Must-not-fire set.
## interpRescue'd results raise ONLY the rescue warning (certificate skipped
## on rescued results by construction -- rescued mmA carries no refLocal).

## Basin-bistable (F4 second manifestation, M4/MP2): the container solve
## goes singular -> interpRescue -> rescue warning; HOMEDESKTOP converges
## directly at the floor -> near-eps warning. Both are honest; either text
## passes, but SOME machine-precision warning must fire.
w_x2 <- tryCatch(minimaxApprox(function(x) x^2, -1, 1, 2L),
                 warning = function(w) conditionMessage(w))
expect_true(grepl("NOT technically a Remez result", w_x2, fixed = TRUE) ||
              grepl("very near machine double precision", w_x2,
                    fixed = TRUE))

## Floor-regime converged fit (exp on [5, 6], deg 10, relaxed convrat; M6
## demonstration case): the 1.05 grid ratio there is linear-solve noise below
## the floor gate, so the certificate must stay silent (asserted where the
## case lives, test_Chebyshev.R; here assert the driver-level skip directly).
o56 <- minimaxApprox(exp, 5, 6, 10L, basis = "c", opts = list(convrat = 1.03))
expect_false(o56$Warning)
## Runge degree-10 (issue #2 restart path, DP-4b recompute): warning-free
## with the pinned ExpErr.
runge <- function(x) 1 / (1 + (5 * x) ^ 2)
o_runge <- sM(minimaxApprox(runge, -1, 1, 10L, basis = "m"))
expect_false(o_runge$Warning)
expect_equal(o_runge$ExpErr, 0.06592293, tolerance = 1e-7)
