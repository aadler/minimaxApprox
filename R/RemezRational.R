# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix
remRatMat <- function(x, E, y, nD, dD) {
  n <- length(x)
  altSgn <- (-1) ^ (seq_len(n) - 1)
  Evctr <- altSgn * E
  yvctr <- -(y + Evctr)
  AMat <- vanderMat(x, nD)
  BMat <- vanderMat(x, dD)[, -1] * yvctr
  cbind(AMat, BMat, -altSgn, deparse.level = 0)
}

# Function to calculate coefficients given matrix and known values
remRatCoeffs <- function(x, E, fn, nD, dD) {
  y <- callFun(fn, x)
  P <- solve(remRatMat(x, E, y, nD, dD), y)
  RR <- list(a = P[seq_len(nD + 1)],            # Works even if nD = 0
             b = c(1, P[seq_len(dD) + nD + 1]), # Works even if dD = 0
             E = P[length(P)])

  if (sum(lengths(RR)) != (length(x) + 1)) { # Not sure this is needed anymore
    stop("Catastrophic Error. Result vector not of required length.")
  }
  RR
}

# Function to calculate value of minimax rational approximation at x given a & b
remRatFunc <- function(x, a, b)  polyCalc(x, a) / polyCalc(x, b)

# Function to calculate error between known and calculated values
remRatErr <- function(x, a, b, fn, absErr) {
  y <- callFun(fn, x)
  (remRatFunc(x, a, b) - y) / if (absErr) 1 else abs(y)
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
remRatRoots <- function(x, a, b, fn, absErr) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remRatErr, interval = intv, a = a, b = b, fn = fn,
                             absErr = absErr),
                     error = function(cond) simpleError(trimws(cond$message)))

    # If there is no root in the interval, take the endpoint closest to zero
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
remRatSwitch <- function(r, l, u, a, b, fn, absErr) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remRatErr(l, a, b, fn, absErr)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remRatErr, interval = intv,
                                a = a, b = b, fn = fn,
                                absErr = absErr, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    if (inherits(extrma, "simpleError")) {
      endPtErr <- remRatErr(intv, a, b, fn, absErr)
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
    E <- remRatErr(p, a, b, fn, absErr)

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
remRat <- function(fn, lower, upper, numerd, denomd, absErr, xi = NULL,
                   opts = list()) {

  # Handle configuration options
  if ("maxiter" %in% names(opts)) {
    maxiter <- opts$maxiter
  } else {
    maxiter <- 100L
  }

  if ("miniter" %in% names(opts)) {
    miniter <- opts$miniter
  } else {
    miniter <- 10L
  }

  if ("conviter" %in% names(opts)) {
    conviter <- opts$conviter
    # If passing conviter then assume maxiter wants at least that much too
    maxiter <- max(maxiter, conviter)
  } else {
    conviter <- 10L
  }

  if ("showProgress" %in% names(opts)) {
    showProgress <- opts$showProgress
  } else {
    showProgress <- FALSE
  }

  if ("cnvgRatio" %in% names(opts)) {
    cnvgRatio <- opts$cnvgRatio
  } else {
    cnvgRatio <- 1 + 1e-11 # Per Cody (1968) can reasonably expect 12 signdig
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- 1e-14
  }

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
  convergeErr <- function(x) {
    E <- 0
    j <- 0L
    repeat {
      if (j > maxiter) break
      RR <- remRatCoeffs(x, E, fn, numerd, denomd)
      if (abs(RR$E - E) <= tol) break
      E <- (RR$E + E) / 2
      j <- j + 1
    }
    RR
  }

  RR <- convergeErr(x)
  errs_last <- remRatErr(x, RR$a, RR$b, fn, absErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    if (i >= maxiter) break
    i <- i + 1L
    r <- remRatRoots(x, RR$a, RR$b, fn, absErr)
    x <- remRatSwitch(r, lower, upper, RR$a, RR$b, fn, absErr)
    RR <- convergeErr(x)
    dngr <- checkDenom(RR$b, lower, upper)
    if (!is.null(dngr)) {
      stop("The ", denomd, " degree polynomial in the denominator has a zero ",
           "at ", fC(dngr), " which makes rational approximation perilous ",
           "over the interval [", fC(lower), ", ", fC(upper), "].")
    }
    errs <- remRatErr(x, RR$a, RR$b, fn, absErr)
    mxae <- max(abs(errs))
    expe <- abs(RR$E)

    if (showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae),
              " Ratio: ", fC(mxae / expe), " Diff:", fC(abs(mxae - expe)))
    }

    # Check for convergence
    if (isConverged(errs, expe, cnvgRatio, tol) && i >= miniter) {
      converged <- TRUE
      break
    }

    # Check that solution is evolving. If solution is not evolving then further
    # iterations will just not help.
    if (all(errs / errs_last <= cnvgRatio) ||
        all(abs(errs - errs_last) <= tol)) {
      unchanging_i <- unchanging_i + 1L
      if (unchanging_i >= conviter) {
        unchanged <- TRUE
        break
      }
    }
    errs_last <- errs
  }

  gotWarning <- FALSE

  if (i >= maxiter && !converged) {
    warning("Convergence to requested ratio and tolerance not acheived in ",
            i, " iterations.\n", "The ratio is ", fC(mxae / expe),
            " times expected and the difference is ", fC(abs(mxae - expe)),
            " from the expected.")
    gotWarning <- TRUE
  }

  if (unchanged && !converged) {
    warning("Convergence to requested ratio and tolerance not acheived in ",
            i, " iterations.\n", unchanging_i, " succesive calculated errors ",
            "were too close to each other to warrant further iterations.\n",
            "The ratio is ", fC(mxae / expe), " times expected and the ",
            "difference is ", fC(abs(mxae - expe)),
            " from the expected.")
    gotWarning <- TRUE
  }

  if (all(abs(errs) < 10 * .Machine$double.eps)) {
    message("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
  }

  ret <- list(a = RR$a, b = RR$b, EE = expe, OE = mxae, iterations = i, x = x,
              Warning = gotWarning)
  attr(ret, "type") <- "Rational"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "absErr") <- absErr
  attr(ret, "tol") <- tol
  attr(ret, "cnvgRatio") <- cnvgRatio
  class(ret) <- c("MiniMaxApprox", class(ret))
  ret
}
