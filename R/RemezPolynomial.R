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
  (polyCalc(x, a) - y) / if (absErr) 1 else y
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
remPolyRoots <- function(x, a, fn, absErr) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remPolyErr, interval = intv, a = a, fn = fn,
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
remPolySwitch <- function(r, l, u, a, fn, absErr) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remPolyErr(l, a, fn, absErr)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remPolyErr, interval = intv, a = a, fn = fn,
                                absErr = absErr, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    # Check if no root in interval and return appropriate endpoint
    if (inherits(extrma, "simpleError")) {
      endPtErr <- remPolyErr(intv, a, fn, absErr)
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
    E <- remPolyErr(p, a, fn, absErr)

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

  if ("convRatio" %in% names(opts)) {
    convRatio <- opts$convRatio
  } else {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    convRatio <- 1.000000001
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
  errs_last <- remPolyErr(x, PP$a, fn, absErr)
  converged <- FALSE
  unchanged <- FALSE
  unchanging_i <- 0L
  i <- 0L
  repeat {
    # Check for maxiter
    if (i >= maxiter) break
    i <- i + 1L
    r <- remPolyRoots(x, PP$a, fn, absErr)
    x <- remPolySwitch(r, lower, upper, PP$a, fn, absErr)
    PP <- remPolyCoeffs(x, fn)
    errs <- remPolyErr(x, PP$a, fn, absErr)
    mxae <- max(abs(errs))
    expe <- abs(PP$E)

    if (showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae))
    }

    # Check for convergence
    if (isConverged(errs, expe, convRatio, tol) && i >= miniter) {
      converged <- TRUE
      break
    }

    # Check that solution is evolving. If solution is not evolving then further
    # iterations will just not help.
    if (all(errs / errs_last <= convRatio) ||
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

  if (mxae < 10 * .Machine$double.eps) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
    gotWarning <- TRUE
  }

  ret <- list(a = PP$a, EE = mxae, OE = mxae,  iterations = i, x = x,
              Warning = gotWarning)
  attr(ret, "type") <- "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "absErr") <- absErr
  attr(ret, "tol") <- tol
  attr(ret, "convRatio") <- convRatio
  class(ret) <- c("minimaxApprox", class(ret))

  ret
}
