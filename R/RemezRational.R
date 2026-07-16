# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# F13 (tolerance half): rational QR rank-detection tolerance. Kept DISTINCT
# from the polynomial path's QRTOLPOLY (1e-14, defined in RemezPolynomial.R)
# on purpose -- see the full rationale there. The rational path has no
# degree-restart machinery, so its rank threshold is not load-bearing in the
# same way, but the two are documented as an intentional pair rather than
# silently divergent magic numbers.
QRTOLRAT <- .Machine$double.eps

# Function to create augmented Vandermonde or Chebyshev matrix for rational
# approximation.
ratMat <- function(x, E, y, nD, dD, relErr, basis, l, u) {
  altSgn <- (-1) ^ (seq_along(x) - 1L)
  # For relative error, need to weight the E by f(x).
  if (relErr) altSgn <- altSgn * y
  altE <- altSgn * E
  yvctr <- -(y + altE)
  # M6 (F5 Option A): l/u required (see evalFunc in shared.R for rationale).
  z <- if (basis == "c") chebMap(x, l, u) else x
  aMat <- if (basis == "m") vanderMat(x, nD) else chebMat(z, nD)
  bMat <- (if (basis == "m") vanderMat(x, dD) else chebMat(z, dD))[, -1L] *
    yvctr
  cbind(aMat, bMat, -altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values.
ratCoeffs <- function(x, E, fn, nD, dD, relErr, basis, l, u, zt) {
  y <- callFun(fn, x)
  P <- ratMat(x, E, y, nD, dD, relErr, basis, l, u)
  PP <- tryCatch(solve(P, y),
                 error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y,
                                                  tol = QRTOLRAT)
  list(a = checkIrrelevant(PP[seq_len(nD + 1L)], l, u, zt, basis),
       b = checkIrrelevant(c(1, PP[seq_len(dD) + nD + 1L]), l, u, zt, basis),
       E = PP[length(PP)])
}

# Main function to calculate and return the minimax rational approximation.
remRat <- function(fn, lower, upper, numerd, denomd, relErr, basis, xi, opts) {

  # Set ZeroBasis relErr flag
  relErrZeroBasis <- FALSE

  # Initial x's
  nodeCount <- numerd + denomd + 2L
  if (is.null(xi)) {
    x <- chebNodes(nodeCount, lower, upper)
  } else if (length(xi) == nodeCount) {
    x <- xi
  } else {
    stop("Given the requested degrees for numerator and denominator, the ",
         "x-vector needs to have ", nodeCount, " elements.")
  }

  # Since E is initially a guess, we need to iterate solving the system of
  # equations until E converges. This function remains *inside* of remRat;
  # therefore, Everything but "x" is previously defined and constant inside the
  # main remRat function and does not need to be passed.
  convergeErr <- function(x) {
    E <- 0
    j <- 0L
    repeat {
      if (j >= opts$maxiter) break
      j <- j + 1L
      RR <- ratCoeffs(x, E, fn, numerd, denomd, relErr, basis, lower, upper,
                      opts$ztol)
      if (abs(RR$E - E) <= opts$tol) break
      E <- (RR$E + E) / 2
    }

    RR
  }

  RR <- convergeErr(x)
  errs_last <- remErr(x, RR, fn, relErr, basis, lower, upper)
  converged <- unchanged <- FALSE
  unchanging_i <- i <- 0L
  repeat {
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, RR, fn, relErr, basis, lower, upper)
    x <- switchX(r, lower, upper, RR, fn, relErr, basis, x)
    relErrZeroBasis <- relErrZeroBasis || attr(x, "ZeroBasis")
    RR <- convergeErr(x)
    dngr <- checkDenom(RR$b, lower, upper, basis)
    if (!is.null(dngr)) {
      stop("The ", denomd, " degree polynomial in the denominator has a zero ",
           "at ", fC(dngr), " which makes rational approximation perilous ",
           "over the interval [", fC(lower), ", ", fC(upper), "]. Increasing ",
           "the numerator or denominator degree by 1 sometimes allows ",
           "convergence.")
    }
    errs <- remErr(x, RR, fn, relErr, basis, lower, upper)
    mxae <- max(abs(errs))
    expe <- abs(RR$E)

    if (opts$showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae),
              " Ratio: ", fC(mxae / expe), " Diff:", fC(abs(mxae - expe)))
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

  # CP-1: reference-local certificate at the converged exit only. See
  # refLocalCheck in shared.R for mechanism and skip conditions.
  rl <- if (converged) {
    refLocalCheck(RR, fn, relErr, basis, lower, upper, expe)
  } else {
    list(gridSup = NA_real_, refLocal = FALSE)
  }

  list(a = RR$a, b = RR$b, expe = expe, mxae = mxae, i = i, x = x,
       converged = converged, unchanged = unchanged,
       unchanging_i = unchanging_i, zeroBasisError = relErrZeroBasis,
       refLocal = rl$refLocal, gridSup = rl$gridSup)
}
