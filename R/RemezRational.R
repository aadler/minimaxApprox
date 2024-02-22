# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix for rational approximation.
ratMat <- function(x, E, y, nD, dD, relErr, basis) {
  n <- length(x)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x).
  if (relErr)  altSgn <- altSgn * y
  altE <- altSgn * E
  yvctr <- -(y + altE)
  matFunc <- switch(EXPR = basis, "m" = vanderMat, chebMat)
  aMat <- matFunc(x, nD)
  bMat <- matFunc(x, dD)[, -1L] * yvctr
  cbind(aMat, bMat, -altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values.
ratCoeffs <- function(x, E, fn, nD, dD, relErr, basis, l, u, zt) {
  y <- callFun(fn, x)
  P <- ratMat(x, E, y, nD, dD, relErr, basis)
  PP <- tryCatch(solve(P, y),
                 error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y,
                                                  tol = .Machine$double.eps)
  list(a = checkIrrelevant(PP[seq_len(nD + 1L)], l, u, zt),
       b = checkIrrelevant(c(1, PP[seq_len(dD) + nD + 1L]), l, u, zt),
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
  } else {
    x <- xi
    if (length(xi) != nodeCount) {
      stop("Given the requested degrees for numerator and denominator, ",
           "the x-vector needs to have ", nodeCount, " elements.")
    }
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
  errs_last <- remErr(x, RR, fn, relErr, basis)
  converged <- unchanged <- FALSE
  unchanging_i <- i <- 0L
  repeat {
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, RR, fn, relErr, basis)
    x <- switchX(r, lower, upper, RR, fn, relErr, basis)
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
    errs <- remErr(x, RR, fn, relErr, basis)
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

  list(a = RR$a, b = RR$b, expe = expe, mxae = mxae, i = i, x = x,
       converged = converged, unchanged = unchanged,
       unchanging_i = unchanging_i, zeroBasisError = relErrZeroBasis)
}
