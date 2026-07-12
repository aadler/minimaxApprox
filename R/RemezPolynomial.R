# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# F13 (tolerance half): the QR rank-detection tolerance for the polynomial
# least-squares fallback in polyCoeffs (and re-used by the F4 interpolant
# rescue below). This is DELIBERATELY LARGER than the rational path's
# tolerance (QRTOLRAT = .Machine$double.eps, defined in RemezRational.R). The
# two are NOT aligned, on purpose: qr.solve's `tol` is the threshold at which a
# solve is declared rank-deficient, and that decision is what drives the
# degree-restart machinery in minimaxApprox() (a singular polynomial solve is
# the SIGNAL to try degree n + 1). Loosening this toward eps would make
# singularity fire less often and shift when restarts trigger -- the same
# "load-bearing singular solve" hazard class documented in the M3 record. The
# container suite is invariant to aligning them, but that has not been verified
# across BLAS/LAPACK platforms, so the values are kept distinct and named
# rather than merged. If a future master pass verifies invariance on all target
# platforms, the two constants may be unified.
QRTOLPOLY <- 1e-14

# Function to create augmented Vandermonde or Chebyshev matrix for polynomial
# approximation.
polyMat <- function(x, y, relErr, basis) {
  n <- length(x)
  matFunc <- switch(EXPR = basis, m = vanderMat, chebMat)
  A <- matFunc(x, n - 2L)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x).
  if (relErr) altSgn <- altSgn * y
  cbind(A, altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values.
polyCoeffs <- function(x, fn, relErr, basis, l, u, zt) {
  y <- callFun(fn, x)
  P <- polyMat(x, y, relErr, basis)
  PP <- tryCatch(solve(P, y),
                 error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y, tol = QRTOLPOLY)
  list(a = checkIrrelevant(PP[-length(PP)], l, u, zt, basis),
       E = PP[length(PP)])
}

# Main function to calculate and return the minimax polynomial approximation.
remPoly <- function(fn, lower, upper, degree, relErr, basis, opts) {

  # Set ZeroBasis relErr flag
  relErrZeroBasis <- FALSE

  # Initial x's
  x <- chebNodes(degree + 2L, lower, upper)

  # Initial Polynomial Guess
  PP <- polyCoeffs(x, fn, relErr, basis, lower, upper, opts$ztol)
  errs_last <- remErr(x, PP, fn, relErr, basis)
  converged <- unchanged <- FALSE
  unchanging_i <- i <- 0L
  repeat {
    # Check for maxiter
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, PP, fn, relErr, basis)
    x <- switchX(r, lower, upper, PP, fn, relErr, basis)
    relErrZeroBasis <- relErrZeroBasis || attr(x, "ZeroBasis")
    PP <- polyCoeffs(x, fn, relErr, basis, lower, upper, opts$ztol)
    errs <- remErr(x, PP, fn, relErr, basis)
    mxae <- max(abs(errs))
    expe <- abs(PP$E)

    if (opts$showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae))
    }

    # Check for convergence
    if (isConverged(errs, expe, opts$convrat, opts$tol) && i >= opts$miniter) {
      converged <- TRUE
      break
    }

    # Check that solution is evolving. If solution is not evolving then further
    # iterations will not help.
    if (isUnchanging(errs, errs_last, opts$convrat, opts$tol)) {
      unchanging_i <- unchanging_i + 1L
      if (unchanging_i >= opts$conviter) {
        unchanged <- TRUE
        break
      }
    }

    errs_last <- errs
  }

  list(a = PP$a, expe = expe, mxae = mxae, i = i, x = x, converged = converged,
       unchanged = unchanged, unchanging_i = unchanging_i,
       zeroBasisError = relErrZeroBasis)
}

# F4 rescue: exactly-representable / machine-precision-resolved interpolant.
#
# WHAT THIS IS. A single plain interpolation, NOT a Remez step. It solves the
# square system p(x_i) = f(x_i) through degree + 1 Chebyshev nodes for the
# degree + 1 coefficients -- no alternating-sign leveled-error (E) column, no
# reference exchange, no iteration.
#
# WHY IT EXISTS. It is reached ONLY after both the degree-n and degree-(n+1)
# Remez solves in minimaxApprox() have failed singular. That double failure
# has two possible causes: (a) the requested function is (numerically) exactly
# a polynomial of degree <= n, so the augmented Remez system has no unique
# solution because there is no minimax problem left to solve (the answer is f
# itself, leveled error 0); or (b) a genuine non-convergence unrelated to
# representability. This routine distinguishes the two by MEASUREMENT: it
# builds the interpolant and checks its error on a dense probe grid. In case
# (a) the interpolant reproduces f to the machine-precision floor; in case (b)
# it does not.
#
# WHY THE E COLUMN IS ABSENT AND relErr STILL WORKS. The E column in polyMat()
# is scaled by y precisely to co-solve for the leveled error inside the Remez
# system; it is a device for that extra unknown, not for "computing relative
# error". Here there is no E to solve for. The interpolation target is f(x_i)
# regardless of error mode. relErr enters ONLY in the measurement below
# ((p - f)/f vs p - f), never in the solve.
#
# WHAT IT RETURNS / MINIMAX STATUS. If the probe error is at the machine floor
# it returns the interpolant (caller flags it with a warning that it is not a
# Remez result). For an exactly-representable f this interpolant equals f and
# is the true minimax solution (E = 0). For an f merely resolved to precision
# by the requested degree, it is minimax only to within floating point: the
# gap to the true minimax polynomial is below what a double can represent, so
# no distinct better solution can be exhibited. If the probe error is NOT at
# the floor, the singularity had cause (b): return NULL and let the caller
# raise the original hard error unchanged -- the rescue must never mask a
# genuine non-convergence with a suboptimal interpolant.
#
# THRESHOLD. 10 * eps, matching the existing near-eps warning's convention in
# minimaxApprox(); scaled by max(1, max|f|) on the grid in absolute mode so the
# floor tracks the function's magnitude (a dimensionless 10 * eps in relative
# mode). relErr with a zero of f on the grid makes the relative criterion
# meaningless (division by ~0): do NOT rescue, fall through.
#
# Returns list(a, x, err) on rescue, or NULL to fall through to the hard error.
interpRescue <- function(fn, lower, upper, degree, relErr, basis) {
  # degree + 1 interpolation nodes (the initial Remez reference is degree + 2;
  # here we need an (degree + 1)-point square interpolation system).
  x <- chebNodes(degree + 1L, lower, upper)
  y <- callFun(fn, x)
  matFunc <- switch(EXPR = basis, m = vanderMat, chebMat)
  A <- matFunc(x, degree)

  # Same solve -> qr.solve fallback shape as polyCoeffs, but on the UNaugmented
  # interpolation matrix, which is generally far better conditioned than the
  # Remez augmented matrix that just failed. QRTOLPOLY keeps the rank threshold
  # consistent with the polynomial path (F13).
  a <- tryCatch(solve(A, y),
                error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(a, "simpleError")) {
    a <- tryCatch(qr.solve(A, y, tol = QRTOLPOLY),
                  error = function(cond) simpleError(trimws(cond$message)))
  }
  # If even the plain interpolation is singular, there is nothing to rescue.
  if (inherits(a, "simpleError")) return(NULL)

  # Measure the interpolant against f on a dense probe grid. Grid density was
  # verified stable (probe error flat from 1e3 to 5e4 points on all F4 cases).
  i_grid <- seq(lower, upper, length.out = 2001L)
  fg <- callFun(fn, i_grid)
  calcFn <- switch(EXPR = basis, m = polyCalc, chebCalc)
  pg <- calcFn(i_grid, a)

  if (relErr) {
    # Relative criterion is undefined where f == 0: do not rescue.
    if (any(fg == 0)) return(NULL)
    err <- max(abs((pg - fg) / fg))
    thresh <- 10 * .Machine$double.eps
  } else {
    err <- max(abs(pg - fg))
    thresh <- 10 * .Machine$double.eps * max(1, max(abs(fg)))
  }

  # Not at the machine floor: the singularity had a different cause. Fall
  # through so the caller raises the original "neither converged" error.
  if (err > thresh) return(NULL)

  list(a = a, x = x, err = err)
}

