# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix for rational approximation
ratMat <- function(x, E, y, nD, dD, relErr) {
  n <- length(x)
  altSgn <- (-1) ^ (seq_len(n) - 1L)
  # For relative error, need to weight the E by f(x)
  if (relErr)  altSgn <- altSgn * y
  altE <- altSgn * E
  yvctr <- -(y + altE)
  aMat <- vanderMat(x, nD)
  bMat <- vanderMat(x, dD)[, -1L] * yvctr
  cbind(aMat, bMat, -altSgn, deparse.level = 0L)
}

# Function to calculate coefficients given matrix and known values
ratCoeffs <- function(x, E, fn, nD, dD, relErr) {
  y <- callFun(fn, x)
  P <- ratMat(x, E, y, nD, dD, relErr)
  PP <- tryCatch(solve(P, y),
                 error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(PP, "simpleError")) PP <- qr.solve(P, y, tol = 1e-12)
  list(a = PP[seq_len(nD + 1L)],            # Works even if nD = 0
       b = c(1, PP[seq_len(dD) + nD + 1L]), # Works even if dD = 0
       E = PP[length(PP)])
}

# Main function to calculate and return the minimax rational approximation
remRat <- function(fn, lower, upper, numerd, denomd, relErr, xi, opts) {

  # Initial x's
  if (is.null(xi)) {
    x <- chebNodes(numerd + denomd + 2L, lower, upper)
  } else {
    x <- xi
    if (length(xi) != numerd + denomd + 2L) {
      stop("Given the requested degrees for numerator and denominator, ",
           "the x-vector needs to have ", numerd + denomd + 2, " elements.")
    }
  }

  # Since E is initially a guess we need to iterate solving the system of
  # equations until E converges. Function may remain inside of remRat.
  # Everything but "x" is previously defined and constant inside the main remRat
  # function and thus does not need to be passed.
  convergeErr <- function(x) {
    E <- 0
    j <- 0L
    repeat {
      if (j >= opts$maxiter) break
      RR <- ratCoeffs(x, E, fn, numerd, denomd, relErr)
      if (abs(RR$E - E) <= opts$tol) break
      E <- (RR$E + E) / 2
      j <- j + 1
    }
    RR
  }

  RR <- convergeErr(x)
  errs_last <- remErr(x, RR, fn, relErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    if (i >= opts$maxiter) break
    i <- i + 1L
    r <- findRoots(x, RR, fn, relErr)
    x <- switchX(r, lower, upper, RR, fn, relErr)
    RR <- convergeErr(x)
    dngr <- checkDenom(RR$b, lower, upper)
    if (!is.null(dngr)) {
      stop("The ", denomd, " degree polynomial in the denominator has a zero ",
           "at ", fC(dngr), " which makes rational approximation perilous ",
           "over the interval [", fC(lower), ", ", fC(upper), "]. Increasing ",
           "the numerator or denominator degree by 1 sometimes allows ",
           "convergence.")
    }
    errs <- remErr(x, RR, fn, relErr)
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
    # iterations will just not help.
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

  list(a = RR$a, b = RR$b, expe = expe, mxae = mxae, i = i, x = x,
       converged = converged, unchanged = unchanged,
       unchanging_i = unchanging_i)
}
