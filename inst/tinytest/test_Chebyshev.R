# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# chebMat
x <- c(-1.5, 0, 1.5)
nx <- length(x)
k <- 3
control <- matrix(c(rep(1, nx), x, 2 * x ^ 2 - 1, 4 * x ^ 3 - 3 * x), ncol = 4L,
                  byrow = FALSE)
expect_equal(minimaxApprox:::chebMat(x, k), control, tolerance = tol)

# chebCalc
a <- 2:5
x <- c(-0.5, 1.5)
control <- c(2 * 1 + 3 * -0.5 + 4 * (2 * (-0.5) ^ 2 - 1) +
               5 * (4 * (-0.5) ^ 3 - 3 * (-0.5)),
             2 * 1 + 3 * 1.5 + 4 * (2 * 1.5 ^ 2 - 1) +
               5 * (4 * 1.5 ^ 3 - 3 * 1.5))
expect_equal(minimaxApprox:::chebCalc(x, a), control, tolerance = tol)

# evalFunc using Chebyshev
## Uses a and x from immediately previous
## M6: evalFunc's l/u are now required. -1,1 is the identity map, so this
## remains exactly the raw-x pass-through check it always was.
## Polynomial
R <- list(a = a)
expect_equal(minimaxApprox:::evalFunc(x, R, "c", -1, 1), control,
             tolerance = tol)

## Rational
b <- c(-0.5, -1, 2)
controlD <- c(-0.5 * 1 - (-0.5) + 2 * (2 * (-0.5) ^ 2 - 1),
              -0.5 * 1 - 1.5 + 2 * (2 * 1.5 ^ 2 - 1))
control <- control / controlD
R <- list(a = a, b = b)
expect_equal(minimaxApprox:::evalFunc(x, R, "c", -1, 1), control,
             tolerance = tol)

# cheb2mon
## Polynomial
A <- minimaxApprox(exp, 0, 1, 4L)
B <- minimaxApprox(exp, 0, 1, 4L, basis = "m")
expect_equal(A$aMono, B$a, tolerance = tol)

## Rational
A <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L))
B <- minimaxApprox(function(x) gamma(x + 1), 1, 2, c(3L, 3L), basis = "m")

expect_equal(A$aMono, B$a, tolerance = tol)
expect_equal(A$bMono, B$b, tolerance = tol)

# --- F1 regression: chebCalc_c large-vector VLA/R_alloc fix -----------------
# Prior to the fix, src/Chebyshev.c's chebCalc_c allocated its working matrix
# as a C stack VLA (`double cMat[m * n]`); a several-million-point evaluation
# overflowed the C call stack and crashed the R session ("segfault from C
# stack overflow"). The fix (R_alloc) must handle this cleanly. One large call
# is used to keep runtime reasonable; a handful of spot-check expectations
# confirm correctness rather than just "didn't crash".
mmC <- minimaxApprox(exp, -1, 1, 5L, basis = "Chebyshev")
mmM <- minimaxApprox(exp, -1, 1, 5L, basis = "monomial")
xBig <- seq(0, 1, length.out = 5e6)
resC <- minimaxEval(xBig, mmC)
resM <- suppressMessages(minimaxEval(xBig, mmM))
expect_length(resC, 5e6)
expect_true(all(is.finite(resC)))
# Chebyshev and monomial representations of the same degree-5 approximation
# must agree to floating tolerance at every point on the grid.
expect_equal(resC, resM, tolerance = 1e-8)
# Spot-check endpoints and midpoint against direct exp() at loose tolerance
# (the polynomial is a degree-5 minimax fit, not exp itself).
expect_equal(resC[c(1, 2.5e6, 5e6)], exp(xBig[c(1, 2.5e6, 5e6)]),
             tolerance = 1e-3)

# --- F12 regression: 2^31-cell guard -----------------------------------
# chebMat_c and chebCalc_c share the same check_cell_cap() guard, so
# exercising it through chebMat_c (cheap: only the scalar degree argument
# need be large, no large vector allocation required) validates the shared
# code path used by both entry points.
expect_error(minimaxApprox:::chebMat(1:2, 2147483647),
             pattern = "exceeds the supported 2^31-cell limit",
             fixed = TRUE)

# --- M6 regression: mapped Chebyshev basis (F5 Option A) --------------------

# chebMap/composeAffine identity on [-1, 1]. This is relied upon throughout
# the package to guarantee every Chebyshev-basis result on [-1, 1] is
# bitwise-unchanged from prior releases.
x <- c(-1, -0.37, 0, 0.6, 1)
expect_identical(minimaxApprox:::chebMap(x, -1, 1), x)
cf <- c(1, -2, 3, -4, 5)
expect_identical(minimaxApprox:::composeAffine(cf, -1, 1), cf)
## Middle coefficient exactly 0, range off [-1,1] (so the identity fast path
## is not taken): exercises the `if (ck == 0) next` skip, otherwise only
## reachable incidentally depending on which coefficients a given Remez fit
## happens to produce. Deterministic, no Remez iteration/BLAS involved.
## z = (2x - (0+2))/(2-0) = x - 1 on [0, 2]; p(z) = 1 + 0*z + 3*z^2, so
## p(x-1) = 1 + 3*(x-1)^2 = 3x^2 - 6x + 4 -> coefficients [4, -6, 3].
expect_identical(minimaxApprox:::composeAffine(c(1, 0, 3), 0, 2), c(4, -6, 3))

# F5 payoff case: exp on [5, 6], deg 10, basis "c". Pre-M6 this converged
# (with a warning) only to ObsErr ~1.8e-10 against unmapped-basis kappa
# ~1.4e24; the true minimax value (measured via an independent barycentric
# prototype, master doc Section 3) is E ~2.9e-12. Tolerance is loose (relative
# to the ~60x improvement being demonstrated) to stay robust across BLAS/
# platform differences, per the project's platform-fragility convention.
r5_6 <- suppressWarnings(minimaxApprox(exp, 5, 6, 10L, basis = "c"))
expect_true(r5_6$ObsErr < 1e-10)
expect_true(r5_6$ObsErr > 1e-13)
# A default-opts Warning is still EXPECTED here (see NEWS/session record): the
# residual oscillation is a floating-point noise floor of the classical
# linear-solve Remez iteration at this magnitude (relative error ~34*eps),
# not a basis-conditioning artifact -- fixing it is out of M6's scope
# (flagged for Master Pass 2). A modestly relaxed convrat (opts-only, NOT a
# package default) reaches the same accuracy without hitting maxiter,
# demonstrating the underlying fit is genuinely converged, not merely
# "improved but still broken".
r5_6_relaxed <- minimaxApprox(exp, 5, 6, 10L, basis = "c",
                              opts = list(convrat = 1.03))
expect_false(r5_6_relaxed$Warning)
expect_equal(r5_6_relaxed$ObsErr, r5_6$ObsErr, tolerance = 0.05)

# Composition accuracy (composeAffine / aMono), normal case: modest degree,
# range close to [-1, 1]. Must stay tight.
rN <- suppressWarnings(minimaxApprox(exp, 0, 1, 4L, basis = "c"))
gridN <- seq(0, 1, length.out = 1001L)
chebN <- minimaxApprox:::chebCalc(minimaxApprox:::chebMap(gridN, 0, 1), rN$a)
monoN <- minimaxApprox:::polyCalc(gridN, rN$aMono)
expect_true(max(abs(chebN - monoN)) < 1e-10)

# Composition accuracy, documented caveat case: high degree, wide range far
# from [-1, 1]. This is the aMono accuracy limit documented in
# man/MiniMaxApprox.Rd (Mapped Chebyshev Basis section) -- the composition
# error here is expected to be non-trivial (and can exceed the underlying
# fit's own accuracy), so this test PINS DOWN and documents that magnitude
# rather than asserting tight accuracy that does not hold for this case.
rW <- suppressWarnings(minimaxApprox(sin, 10, 20, 15L, basis = "c"))
gridW <- seq(10, 20, length.out = 1001L)
chebW <- minimaxApprox:::chebCalc(minimaxApprox:::chebMap(gridW, 10, 20), rW$a)
monoW <- minimaxApprox:::polyCalc(gridW, rW$aMono)
compErrW <- max(abs(chebW - monoW))
expect_true(compErrW > 1e-8)   # confirms the caveat is real, not a fluke
expect_true(compErrW < 1e-4)   # and bounds it, so a future regression here
# (e.g. a broken composeAffine) is still caught

# Cross-basis sanity: exp on [0, 1], degree <= 10, well away from either
# basis's own machine-precision floor (deg 8/10 on this case sit at their own
# floor and were measured to disagree up to ~1% for that reason alone --
# platform-fragile noise, not an M6 regression; excluded here deliberately).
for (d in c(2L, 5L)) {
  rc <- suppressWarnings(minimaxApprox(exp, 0, 1, d, basis = "c"))
  rm <- suppressWarnings(minimaxApprox(exp, 0, 1, d, basis = "m"))
  expect_equal(rc$ObsErr, rm$ObsErr, tolerance = 1e-8)
}

# checkDenom root location stays in raw x under the mapped basis (not the
# internally-mapped z). sin has a root of the denominator polynomial inside
# [0.75*pi, 1.25*pi] at the degree-(2,3) rational fit.
expect_error(minimaxApprox(sin, 0.75 * pi, 1.25 * pi, c(2L, 3L)),
             pattern = "has a zero at 2\\.")
