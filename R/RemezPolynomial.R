# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

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
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y, tol = 1e-14)
  list(a = checkIrrelevant(PP[-length(PP)], l, u, zt), E = PP[length(PP)])
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
