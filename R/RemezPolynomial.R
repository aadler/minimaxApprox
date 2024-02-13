# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix for polynomial approximation.
polyMat <- function(x, y, relErr, monoB) {
  n <- length(x)
  A <- if (monoB) vanderMat(x, n - 2L) else chebMat(x, n - 2L)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x).
  if (relErr) altSgn <- altSgn * y
  cbind(A, altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values.
polyCoeffs <- function(x, fn, relErr, monoB, l, u, zt) {
  y <- callFun(fn, x)
  P <- polyMat(x, y, relErr, monoB)
  PP <- tryCatch(solve(P, y),
                 error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y, tol = 1e-14)
  list(a = checkIrrelevant(PP[-length(PP)], l, u, zt), E = PP[length(PP)])
}

# Main function to calculate and return the minimax polynomial approximation.
remPoly <- function(fn, lower, upper, degree, relErr, monoB, opts) {

  # Set ZeroBasis relErr flag
  relErrZeroBasis <- FALSE

  # Initial x's
  x <- chebNodes(degree + 2L, lower, upper)

  # Initial Polynomial Guess
  PP <- polyCoeffs(x, fn, relErr, monoB, lower, upper, opts$ztol)
  errs_last <- remErr(x, PP, fn, relErr, monoB)
  converged <- unchanged <- FALSE
  unchanging_i <- i <- 0L
  repeat {
    # Check for maxiter
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, PP, fn, relErr, monoB)
    x <- switchX(r, lower, upper, PP, fn, relErr, monoB)
    relErrZeroBasis <- relErrZeroBasis || attr(x, "ZeroBasis")
    PP <- polyCoeffs(x, fn, relErr, monoB, lower, upper, opts$ztol)
    errs <- remErr(x, PP, fn, relErr, monoB)
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
    if (all(errs / errs_last <= opts$convrat) ||
          all(abs(errs - errs_last) <= opts$tol)) {
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
