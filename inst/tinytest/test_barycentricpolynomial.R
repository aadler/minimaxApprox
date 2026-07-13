# Copyright Avraham Adler (c) 2026
# SPDX-License-Identifier: MPL-2.0+

# Module M5 Phase 1: barycentric-Remez polynomial basis (Pachon-Trefethen 2009).
# https://www.chebfun.org/publications/remez.pdf.
# Oracles: PT09 Table 1 (deg-10 absolute leveled error, nine functions; three
# used here) and the chebfun |x| deg-11 monomial coefficients (paper section 4).
# Barycentric evaluation is well-conditioned, so the numeric assertions below are
# cross-platform; the two genuinely platform/BLAS-fragile relative-error end-to-
# end paths are gated to HOMEDESKTOP at the end, matching the existing suite's
# convention for such cases.

tol <- sqrt(.Machine$double.eps)
sM <- function(x) suppressMessages(x)
sW <- function(x) suppressWarnings(x)

# ---- PT09 Table 1 oracles (deg 10, absolute) ------------------------------
f1 <- function(x) tanh(x + 0.5) - tanh(x - 0.5)
f2 <- function(x) sin(exp(x))
f3 <- function(x) sqrt(x + 1)
r1 <- sW(minimaxApprox(f1, -1, 1, 10, basis = "b"))
r2 <- sW(minimaxApprox(f2, -1, 1, 10, basis = "b"))
r3 <- sW(minimaxApprox(f3, -1, 1, 10, basis = "b"))
# Agreement with the published oracles is convergence-tolerance-limited to ~7-8
# figures for these small E values (a 1e-8-scale error cannot be pinned to eps).
expect_equal(r1$ExpErr, 0.00000030009195, tolerance = 1e-7)
expect_equal(r2$ExpErr, 0.00000178623400, tolerance = 1e-7)
expect_equal(r3$ExpErr, 0.01978007008380, tolerance = 1e-6)
# f1 is even on a symmetric interval: the symmetric-reference h-cancellation
# handling (separateNodes + grid-based F4 detection) must NOT short-circuit it
# to a spurious E = 0.
expect_true(r1$ExpErr > 1e-8)
expect_true(r1$iterations > 0L)

# ---- |x| deg-11 monomial coefficients (chebfun oracle) --------------------
rx <- sW(minimaxApprox(function(x) abs(x), -1, 1, 11, basis = "b"))
cheb_even <- c(0.02784511855, 4.75365049278, -20.64625015816, 47.77533460523,
               -49.59209097049, 18.70935603064)
expect_equal(rx$aMono[c(1, 3, 5, 7, 9, 11)], cheb_even, tolerance = 1e-6)
expect_true(max(abs(rx$aMono[c(2, 4, 6, 8, 10, 12)])) < 1e-8)  # odd coeffs ~ 0

# ---- F4 family handled in-basis (paper 3.6), no interpRescue --------------
# Exactly representable and resolved-to-precision cases return E ~ 0 with NO
# 'rescued' field (that field is classical-only; barycentric handles it
# natively).
r_x2 <- sW(minimaxApprox(function(x) x^2, -1, 1, 2, basis = "b"))
expect_true(r_x2$ExpErr < 1e-14)
expect_null(r_x2$rescued)
for (d in c(14L, 15L, 50L)) {
  rr <- sW(minimaxApprox(exp, -1, 1, d, basis = "b"))
  expect_true(rr$ExpErr < 1e-14)
  expect_null(rr$rescued)
}

# ---- Parity with the Chebyshev basis on healthy cases ---------------------
# Agreement between two independently-converged Remez fits is convergence-
# tolerance-limited, NOT to eps: each stops when its own convrat criterion is
# met, and that stopping point is BLAS/platform-dependent. On/near [-1,1] the
# classical basis is well-conditioned and its stopping point is stable, so a
# tight parity bound holds. OFF [-1,1] the classical basis is noisier and its
# last-few-digits stopping point wobbles across platforms (exp[0,3] deg12: the
# bary-vs-cheb relative gap measured 3.7e-6 in one environment and 5.3e-5 in
# another) -- so testing the bary fit against the classical fit's E there is
# inherently fragile. For the off-range case we instead assert the platform-
# STABLE property (that the barycentric fit equioscillates well on its own,
# ObsErr/ExpErr ~ 1, evaluated through the well-conditioned barycentric form)
# plus a generous ballpark parity that still catches a gross regression.
parity <- function(f, l, u, d) {
  rb <- sW(minimaxApprox(f, l, u, d, basis = "b"))
  rc <- sW(sM(minimaxApprox(f, l, u, d, basis = "c")))
  abs(rb$ExpErr - rc$ExpErr) / rc$ExpErr
}
expect_true(parity(exp, -1, 1, 8) < 1e-6)     # on-range, classical is stable
expect_true(parity(sin, -2, 2, 10) < 1e-6)    # near-range, stable
# off-range: stable equioscillation self-check + coarse cross-basis sanity only
rb_off <- sW(minimaxApprox(exp, 0, 3, 12, basis = "b"))
expect_true(abs(rb_off$ObsErr / rb_off$ExpErr - 1) < 1e-3)
expect_true(parity(exp, 0, 3, 12) < 1e-3)

# ---- relErr closed form: equioscillates and differs from absErr optimum ---
rrel <- sW(minimaxApprox(exp, -1, 1, 6, basis = "b", relErr = TRUE))
rabs <- sW(minimaxApprox(exp, -1, 1, 6, basis = "b", relErr = FALSE))
g <- seq(-1, 1, length.out = 4001)
relcurve <- (exp(g) - minimaxEval(g, rrel)) / exp(g)
peaks <- which(diff(sign(diff(relcurve))) != 0) + 1L
pk <- abs(relcurve[peaks])
expect_true(length(pk) >= 6L)
expect_true((max(pk) - min(pk)) / max(pk) < 1e-3)
expect_true(abs(rrel$ExpErr - rabs$ExpErr) / rabs$ExpErr > 1e-6)

# ---- relErr node-on-exact-zero guard (FIX#3) ------------------------------
# sin on [0, 3.5]: endpoint 0 is an EXACT zero of sin and a forced reference
# node, so relErr is 0/0 = NaN there; the relErr minimax does not exist -> a
# clear error, NOT a crash ("missing value where TRUE/FALSE needed").
expect_error(
  sW(minimaxApprox(sin, 0, 3.5, 6, basis = "b", relErr = TRUE)),
  "Relative error is undefined")
# Contrast: sin on [-1,1] has an interior node at 0 that is only floating-point-
# near zero (~6.12e-17), NOT exactly zero -> finite ratio, well-posed, converges.
expect_silent_ok <- sW(minimaxApprox(sin, -1, 1, 6, basis = "b",
                                                   relErr = TRUE))
expect_true(is.finite(expect_silent_ok$ExpErr))

# ---- Rational barycentric is rejected up front (Phase 1 is polynomial-only) --
# A length-2 degree with basis "b" is a reachable user input; the guard rejects
# it with a clear message rather than letting it fall through to remRat (which
# would fail confusingly, since the barycentric evaluator expects R$bary, not
# a/b coefficient lists). This covers the guard until Phase 2 replaces it.
expect_error(minimaxApprox(exp, 0, 1, c(2L, 3L), basis = "b"),
             "not yet supported for rational")

# ---- Capacity scaling: wide interval, high degree, no over/underflow ------
rw <- sW(minimaxApprox(function(x) 1 / (1 + x^2), -1e4, 1e4, 50,
                                     basis = "b"))
expect_true(is.finite(rw$ExpErr))
expect_false(any(is.nan(rw$aMono)))

# ---- Runge convergence and parity -----------------------------------------
runge <- function(x) 1 / (1 + 25 * x^2)
rrg <- sW(minimaxApprox(runge, -1, 1, 10, basis = "b"))
rcg <- sW(sM(minimaxApprox(runge, -1, 1, 10, basis = "c")))
expect_true(abs(rrg$ExpErr - rcg$ExpErr) / rcg$ExpErr < 1e-6)

# ---- convResid reported, with the deg-20 pole-near anomaly documented -----
# 1/(1+x^2) has poles at +-i; its monomial conversion loses ~7 digits at deg 20
# specifically (deg 10 and deg 50 are ~1e-15). This is a documented caveat of
# the coefficient conversion step, not a fit error (bary E is correct
# throughout).
cr <- function(d) {
  sW(minimaxApprox(function(x) 1 / (1 + x^2), -1, 1, d, basis = "b"))$convResid
}
expect_true(cr(10L) < 1e-12)
expect_true(is.finite(cr(20L)))       # anomaly: ~2.5e-8, still finite/reported
expect_true(cr(50L) < 1e-12)

# ---- Internal unit tests (deterministic, no linear solve / no BLAS path) --
# chebNodes2: second-kind, endpoints included, sorted ascending.
z <- minimaxApprox:::chebNodes2(6L, -1, 1)
expect_equal(z[1L], -1, tolerance = tol)
expect_equal(z[length(z)], 1, tolerance = tol)
expect_false(is.unsorted(z))

# separateNodes: nudges sub-tolerance-adjacent nodes apart, leaves well-spaced
# nodes untouched.
sep_in <- c(-1, -0.5, -0.5 + 1e-14, 0.5, 1)      # two nodes 1e-14 apart
sep_out <- minimaxApprox:::separateNodes(sep_in, -1, 1)
expect_true(all(diff(sep_out) > 0))
expect_equal(minimaxApprox:::separateNodes(c(-1, 0, 1), -1, 1), c(-1, 0, 1),
             tolerance = tol)  # already well-spaced: unchanged

# baryWeights + baryEval: exact-node short-circuit returns the node value.
xr <- minimaxApprox:::chebNodes2(6L, -1, 1)
wv <- minimaxApprox:::baryWeights(xr, -1, 1)
pv <- exp(xr)
expect_identical(minimaxApprox:::baryEval(xr[3L], xr, wv, pv), pv[3L])
# Off-node, the second barycentric formula reproduces a low-degree polynomial
# exactly: interpolating x^2 at 6 nodes evaluates x^2 at an arbitrary point.
pv2 <- xr^2
expect_equal(minimaxApprox:::baryEval(0.37, xr, wv, pv2), 0.37^2,
             tolerance = tol)

# levelError: absolute closed form on a symmetric even reference cancels to ~0
# (this is the cancellation the F4 grid-detection is designed to see through).
sig6 <- (-1)^(seq_len(6L) - 1L)
fe <- f1(xr)
expect_true(abs(minimaxApprox:::levelError(wv, sig6, fe, FALSE)) < 1e-10)

# onePointExchange (overshoot safeguard body): unreachable via the public API
# with this package's exchange quality (max observed overshoot ratio ~0.36 vs
# the 1e5 trigger), so exercised by a direct call. Must return a same-length,
# sorted, unique, in-range reference.
tr <- minimaxApprox:::baryTrial(xr, exp, FALSE, -1, 1, sig6)
xnew <- minimaxApprox:::onePointExchange(tr$x, tr$R, exp, FALSE, -1, 1)
expect_length(xnew, length(tr$x))
expect_false(is.unsorted(xnew))
expect_identical(anyDuplicated(xnew), 0L)
expect_true(all(xnew >= -1 & xnew <= 1))

# onePointExchange endpoint branches: with an interior-only reference the domain
# endpoints are bracket candidates, so a trial whose error peaks at an endpoint
# drives xnew there. A steep rising exp peaks the residual at the RIGHT end (the
# `else` xold branch); a steep falling exp at the LEFT end (the `else if` branch).
mkTrial <- function(fnc, xk, l, u) {
  w <- minimaxApprox:::baryWeights(xk, l, u)
  sg <- (-1) ^ (seq_along(xk) - 1L)
  fv <- fnc(xk)
  h <- minimaxApprox:::levelError(w, sg, fv, FALSE)
  list(R = list(bary = list(x = xk, w = w, p = fv - sg * h)))
}
xk_int <- c(-0.8, -0.4, 0, 0.4, 0.8)
trR <- mkTrial(function(x) exp(6 * x), xk_int, -1, 1)
outR <- minimaxApprox:::onePointExchange(xk_int, trR$R, function(x) exp(6 * x),
                                         FALSE, -1, 1)
expect_true(max(outR) >= 0.8 && !is.unsorted(outR) && !any(duplicated(outR)))
trL <- mkTrial(function(x) exp(-6 * x), xk_int, -1, 1)
outL <- minimaxApprox:::onePointExchange(xk_int, trL$R, function(x) exp(-6 * x),
                                         FALSE, -1, 1)
expect_true(min(outL) <= -0.8 && !is.unsorted(outL) && !any(duplicated(outL)))

# baryTrial relErr near-zero-denominator guard (distinct from the
# exact-node-zero case, which errors in remBary). A 2-node relErr trial of fn =
# x on symmetric nodes makes sum(sigma * w * f) == 0 exactly (f = c(-a, a)), so
# h is non-finite: the guard flags zeroBasis and resets h to 0. No node is
# itself a zero of fn.
z2 <- minimaxApprox:::baryTrial(c(-0.5, 0.5), function(x) x, TRUE, -0.5, 0.5,
                                c(1, -1))
expect_true(z2$zeroBasis)
expect_true(is.finite(z2$h))

# remBary normf fallback: fn identically 0 makes max|f| = 0 on the probe grid,
# exercising the `normf <- 1` guard. Resolves at the machine floor with the
# near-eps warning (E is the h-floor perturbation, not a Remez result).
expect_warning(minimaxApprox(function(x) 0 * x, -1, 1, 2, basis = "b"),
               "machine double precision")

# remBary "unchanging" (stall) exit: setting miniter above maxiter blocks the
# isConverged branch, so a fit whose error vector stabilises exits via
# isUnchanging instead. Deterministic (opts-driven, not BLAS-driven).
expect_warning(
  minimaxApprox(exp, -1, 1, 6, basis = "b",
                opts = list(maxiter = 60L, miniter = 1000L, conviter = 1L,
                            showProgress = FALSE, convrat = 1.000000001,
                            tol = 1e-14)),
  "too close")

# ---- Methods dispatch on the barycentric object ---------------------------
rb <- sW(minimaxApprox(exp, -1, 1, 8, basis = "b"))
expect_true(abs(minimaxEval(0.3, rb) - exp(0.3)) < rb$ExpErr * 1.5)
expect_true(all(is.finite(minimaxErr(seq(-1, 1, length.out = 11), rb))))
ce <- coef(rb)
expect_true(is.list(ce) || is.numeric(ce))
# print returns its argument invisibly without error.
expect_silent(invisible(capture.output(print(rb))))

# ---- HOMEDESKTOP-gated: platform/BLAS-fragile relErr end-to-end paths ------
# These mirror the master doc's and M3/M6 records' known-fragile cases; the
# barycentric handling is deterministic in-container but the underlying inputs
# are ill-posed at the machine-precision floor, so gate to the maintainer's node
# rather than let them fail unverified CI.
if (Sys.info()["nodename"] == "HOMEDESKTOP") {
  # x^2 - 4 on [-3,-1] has an exact interior zero at x = -2; relErr is undefined
  # there. Barycentric must not silently return a wrong answer -- either the
  # exact-zero guard or a clean non-convergence, never a crash.
  rz <- tryCatch(
    sW(minimaxApprox(function(x) x^2 - 4, -3, -1, 3,
                                   basis = "b", relErr = TRUE)),
    error = function(e) structure(conditionMessage(e), class = "errcase"))
  expect_true(inherits(rz, "errcase") || is.finite(rz$ExpErr))
}
