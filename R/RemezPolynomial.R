# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix
polyMat <- function(x, y, absErr) {
  n <- length(x)
  A <- vanderMat(x, n - 2L)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x)
  if (!absErr) {
    altSgn <- altSgn * y
  }
  cbind(A, altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values
polyCoeffs <- function(x, fn, absErr) {
  y <- callFun(fn, x)
  PP <- solve(polyMat(x, y, absErr), y)
  list(a = PP[-length(PP)], E = PP[length(PP)])
}

# Main function to calculate and return the minimax polynomial approximation
remPoly <- function(fn, lower, upper, degree, absErr, opts) {

  # Initial x's
  x <- chebNodes(degree + 2L, lower, upper)

  # Initial Polynomial Guess
  PP <- polyCoeffs(x, fn, absErr)
  errs_last <- remErr(x, PP, fn, absErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    # Check for maxiter
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, PP, fn, absErr)
    x <- switchX(r, lower, upper, PP, fn, absErr)
    PP <- polyCoeffs(x, fn, absErr)
    errs <- remErr(x, PP, fn, absErr)
    mxae <- max(abs(errs))
    expe <- abs(PP$E)

    if (opts$showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae))
    }

    # Check for convergence
    if (isConverged(errs, expe, opts$convRatio, opts$tol) &&
        i >= opts$miniter) {
      converged <- TRUE
      break
    }

    # Check that solution is evolving. If solution is not evolving then further
    # iterations will just not help.
    if (all(errs / errs_last <= opts$convRatio) ||
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
       unchanged = unchanged, unchanging_i = unchanging_i)
}
