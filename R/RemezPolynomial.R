# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Function to create augmented Vandermonde matrix
remPolyMat <- function(x) {
  n <- length(x)
  A <- vanderMat(x, n - 2)
  cbind(A, (-1) ^ (seq_len(n) - 1))
}

# Function to calculate coefficients given matrix and known values
remPolyCoeffs <- function(x, fn) {
  PP <- solve(remPolyMat(x), callFun(fn, x))
  list(a = PP[-length(PP)], E = PP[length(PP)])
}

# Function to calculate error between known and calculated values
remPolyErr <- function(x, a, fn, absErr) {
  y <-  callFun(fn, x)
  (polyCalc(x, a) - y) / if (absErr) 1 else abs(y)
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
remPolyRoots <- function(x, a, fn, absErr) {
  # if (all(abs(remPolyErr(x, a, fn, absErr)) <= 5 * .Machine$double.eps)) {
  #   stop("This code only functions to machine double precision. All error ",
  #        "values are too near machine double precision. Please try again ",
  #        "using a lesser degree.")
  # }
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remPolyErr, interval = intv, a = a, fn = fn,
                             absErr = absErr),
                     error = function(cond) simpleError(trimws(cond$message)))
    # If there is no root in the interval, take the lower endpoint
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
remPolySwitch <- function(r, l, u, a, fn, absErr) {
  nodes <- data.frame(lower = c(l, r), upper = c(r, u))
  x <- double(dim(nodes)[1L])
  maximize <- sign(remPolyErr(l, a, fn, absErr)) == 1
  for (i in seq_along(x)) {
    x[i] <- optimize(remPolyErr,
                     lower = nodes$lower[i],
                     upper = nodes$upper[i],
                     a = a, fn = fn, absErr = absErr,
                     maximum = maximize)[[1L]]
    # Test endpoints for max/min
    p <- c(nodes$lower[i], x[i], nodes$upper[i])
    testDF <- data.frame(p = p, E = remPolyErr(p, a, fn, absErr))
    if (maximize) {
      x[i] <- testDF$p[which.max(testDF$E)]
    } else {
      x[i] <- testDF$p[which.min(testDF$E)]
    }

    # Flip maximize
    maximize <- !maximize
  }
  x
}

# Main function to calculate and return the minimax polynomial approximation
remPoly <- function(fn, lower, upper, degree, absErr, opts = list()) {

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

  if ("showProgress" %in% names(opts)) {
    showProgress <- opts$showProgress
  } else {
    showProgress <- FALSE
  }

  if ("cnvgRatio" %in% names(opts)) {
    cnvgRatio <- opts$cnvgRatio
  } else {
    cnvgRatio <- 1.001
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- 1e-14
  }

  # Initial x's
  x <- chebNodes(degree + 2L, lower, upper)

  # Initial Polynomial Guess
  PP <- remPolyCoeffs(x, fn)
  i <- 0L
  repeat {
    i <- i + 1L
    r <- remPolyRoots(x, PP$a, fn, absErr)
    x <- remPolySwitch(r, lower, upper, PP$a, fn, absErr)
    PP <- remPolyCoeffs(x, fn)
    errs <- remPolyErr(x, PP$a, fn, absErr)
    mxae <- max(abs(errs))

    if (showProgress) {
      message("i: ", i, " E: ", fC(PP$E), " maxErr: ", fC(mxae))
    }

    if ((isConverged(errs, abs(PP$E), cnvgRatio, tol) && i >= miniter) ||
        i > maxiter) break
  }

  gotWarning <- FALSE

  if (i >= maxiter) {
    warning("Convergence not acheived in ", maxiter, " iterations.\n",
            "Maximum observed error is ", fC(mxae / PP$E), " times expected.")
    gotWarning <- TRUE
  }

  if (all(abs(errs) < 10 * .Machine$double.eps)) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
    gotWarning <- TRUE
  }

  ret <- list(a = PP$a, EE = abs(PP$E), OE = mxae,  iterations = i, x = x,
              Warning = gotWarning)
  attr(ret, "type") <- "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "absErr") <- absErr
  attr(ret, "tol") <- tol
  attr(ret, "cnvgRatio") <- cnvgRatio
  class(ret) <- c("MiniMaxApprox", class(ret))

  ret
}
