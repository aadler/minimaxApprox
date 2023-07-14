# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix
ratMat <- function(x, E, y, nD, dD, absErr) {
  n <- length(x)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x)
  if (!absErr) {
    altSgn <- altSgn * y
  }
  altE <- altSgn * E
  yvctr <- -(y + altE)
  aMat <- vanderMat(x, nD)
  bMat <- vanderMat(x, dD)[, -1L] * yvctr
  cbind(aMat, bMat, -altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values
ratCoeffs <- function(x, E, fn, nD, dD, absErr) {
  y <- callFun(fn, x)
  P <- solve(ratMat(x, E, y, nD, dD, absErr), y)
  RR <- list(a = P[seq_len(nD + 1)],            # Works even if nD = 0
             b = c(1, P[seq_len(dD) + nD + 1]), # Works even if dD = 0
             E = P[length(P)])
  RR
}

# Function to identify new x positions. This algorithm uses the multi-switch
# paradigm, not the single switch.
remRatSwitch <- function(r, l, u, RR, fn, absErr) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remErr(l, RR, fn, absErr)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remErr, interval = intv, R = RR, fn = fn,
                                absErr = absErr, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    if (inherits(extrma, "simpleError")) {
      endPtErr <- remErr(intv, RR, fn, absErr)
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
    E <- remErr(p, RR, fn, absErr)

    if (maximize) {
      x[i] <- p[which.max(E)]
    } else {
      x[i] <- p[which.min(E)]
    }

    # Test for 0 value at function if relative error
    if (!absErr) {
      if (callFun(fn, x[i]) == 0) {
        stop("Algorithm is choosing basis point where functional value is ",
             "0. Please approximate using absolute, and not relative error.")
      }
    }

    # Flip maximize
    maximize <- !maximize
  }
  x
}

# Main function to calculate and return the minimax rational approximation
remRat <- function(fn, lower, upper, numerd, denomd, absErr, xi = NULL, opts) {

  # Initial x's
  if (is.null(xi)) {
    x <- chebNodes(numerd + denomd + 2, lower, upper)
  } else {
    x <- xi
    if (length(xi) != numerd + denomd + 2) {
      stop("Given the requested degrees for numerator and denominator, ",
           "the x-vector needs to have ", numerd + denomd + 2, " elements.")
    }
  }

  # Since E is initially a guess we need to iterate solving the system of
  # equations until E converges. Function may remain inside of remRat.
  # Everything but "x" is previously defined and constant inside the main remRat
  # function and thus does not need to be passed.
  convergeErr <- function(x, absErr) {
    E <- 0
    j <- 0L
    repeat {
      if (j > opts$maxiter) break
      RR <- ratCoeffs(x, E, fn, numerd, denomd, absErr)
      if (abs(RR$E - E) <= opts$tol) break
      E <- (RR$E + E) / 2
      j <- j + 1
    }
    RR
  }

  RR <- convergeErr(x, absErr)
  errs_last <- remErr(x, RR, fn, absErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, RR, fn, absErr)
    x <- remRatSwitch(r, lower, upper, RR, fn, absErr)
    RR <- convergeErr(x, absErr)
    dngr <- checkDenom(RR$b, lower, upper)
    if (!is.null(dngr)) {
      stop("The ", denomd, " degree polynomial in the denominator has a zero ",
           "at ", fC(dngr), " which makes rational approximation perilous ",
           "over the interval [", fC(lower), ", ", fC(upper), "].")
    }
    errs <- remErr(x, RR, fn, absErr)
    mxae <- max(abs(errs))
    expe <- abs(RR$E)

    if (opts$showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae),
              " Ratio: ", fC(mxae / expe), " Diff:", fC(abs(mxae - expe)))
    }

    # Check for convergence
    if (isConverged(errs, expe, opts$convRatio, opts$tol) && i >= opts$miniter) {
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

  list(a = RR$a, b = RR$b, expe = expe, mxae = mxae, i = i, x = x,
       converged = converged, unchanged = unchanged,
       unchanging_i = unchanging_i)
}
