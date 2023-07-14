# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix
remPolyMat <- function(x, y, absErr) {
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
remPolyCoeffs <- function(x, fn, absErr) {
  y <- callFun(fn, x)
  PP <- solve(remPolyMat(x, y, absErr), y)
  list(a = PP[-length(PP)], E = PP[length(PP)])
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
remPolyRoots <- function(x, PP, fn, absErr) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remErr, interval = intv, R = PP, fn = fn,
                             absErr = absErr),
                     error = function(cond) simpleError(trimws(cond$message)))
    # If there is no root in the interval, take the endpoint closer to 0 (root)
    if (inherits(root, "simpleError")) {
      r[i] <- intv[which.min(abs(intv))]
    } else {
      r[i] <- root$root
    }
  }
  r
}

# Function to identify new x positions. This algorithm uses the multi-switch
# paradigm, not the single switch.
remPolySwitch <- function(r, l, u, PP, fn, absErr) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remErr(l, PP, fn, absErr)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remErr, interval = intv, R = PP, fn = fn,
                                absErr = absErr, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    # Check if no root in interval and return appropriate endpoint
    if (inherits(extrma, "simpleError")) {
      endPtErr <- remErr(intv, PP, fn, absErr)
      if (maximize) {
        x[i] <- intv[which.max(endPtErr)]
      } else {
        x[i] <- intv[which.min(endPtErr)]
      }
    } else {
      x[i] <- extrma[[1L]]
    }

    # Test endpoints for max/min
    p <- c(bottoms[i], x[i], tops[i])
    E <- remErr(p, PP, fn, absErr)

    if (maximize) {
      x[i] <- p[which.max(E)]
    } else {
      x[i] <- p[which.min(E)]
    }

    # Flip maximize
    maximize <- !maximize
  }
  x
}

# Main function to calculate and return the minimax polynomial approximation
remPoly <- function(fn, lower, upper, degree, absErr, opts) {

  # Initial x's
  x <- chebNodes(degree + 2L, lower, upper)

  # Initial Polynomial Guess
  PP <- remPolyCoeffs(x, fn, absErr)
  errs_last <- remErr(x, PP, fn, absErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    # Check for maxiter
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- remPolyRoots(x, PP, fn, absErr)
    x <- remPolySwitch(r, lower, upper, PP, fn, absErr)
    PP <- remPolyCoeffs(x, fn, absErr)
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
