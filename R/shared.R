# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Printing convenience function
fC <- function(x, d = 6L, f = "g", w = -1L) {
  formatC(x, digits = d, format = f, width = w)
}

# Default Chebyeshev nodes
chebNodes <- function(n, a, b) {
  n <- as.integer(n)
  sort(0.5 * (a + b + (b - a) * cos((2 * seq_len(n) - 1) * pi / (2 * n))))
}

# Create Vandermonde matrix to a polynomial degree n, or n + 1 terms
vanderMat <- function(x, n) {
  n1 <- n + 1L
  matrix(rep(x, each = n1) ^ (seq_len(n1) - 1L), ncol = n1, byrow = TRUE)
}

# Call the function being approximated on the points x
callFun <- function(fn, x) {
  if (!is.function(fn)) stop("Unable to parse function.")
  do.call(match.fun(fn), args = list(x = x))
}

# Check that the values passed are oscillating in sign
isOscil <- function(x) all(abs(diff(sign(x))) == 2)

# Calculate the polynomial approximation. Use in numer & denom for rationals
polyCalc <- function(x, a)  drop(vanderMat(x, length(a) - 1) %*% a)

# Function to calculate value of minimax approximation at x given a & b
evalFunc <- function(x, R) {
  ret <- polyCalc(x, R$a)
  if ("b" %in% names(R)) {
    ret <- ret / polyCalc(x, R$b)
  }
  ret
}

# Function to calculate error between known and calculated values
remErr <- function(x, R, fn, relErr) {
    if (relErr) {
      y <- callFun(fn, x)
      (evalFunc(x, R) - y) / y
    } else {
      evalFunc(x, R) - callFun(fn, x)
    }
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
findRoots <- function(x, R, fn, relErr) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remErr, interval = intv, R = R, fn = fn,
                             relErr = relErr),
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
switchX <- function(r, l, u, R, fn, relErr) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remErr(l, R, fn, relErr)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remErr, interval = intv, R = R, fn = fn,
                                relErr = relErr, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    # If no extremum then the take endpoint with "better" value depending if we
    # are maximizing or minimizing.
    if (inherits(extrma, "simpleError")) {
      endPtErr <- remErr(intv, R, fn, relErr)
      if (maximize) {
        x[i] <- intv[which.max(endPtErr)]
      } else {
        x[i] <- intv[which.min(endPtErr)]
      }
    } else {
      x[i] <- extrma[[1L]]
    }

    # Test endpoints for max/min even if an extremum was found.
    p <- c(bottoms[i], x[i], tops[i])
    E <- remErr(p, R, fn, relErr)

    if (maximize) {
      x[i] <- p[which.max(E)]
    } else {
      x[i] <- p[which.min(E)]
    }

    # Test for 0 value at function if relative error
    if (relErr) {
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

# Check Remez iterations for convergence
isConverged <- function(errs, expe, convRatio, tol) {
  aerrs <- abs(errs)
  mxae <- max(aerrs)
  mnae <- min(aerrs)

  # Check observed errors are close enough to expected by ratio or tolerance
  errDistance <- mxae / expe <= convRatio || abs(mxae - expe) <= tol

  # Check observed errors are close enough to each other by ratio or tolerance
  errMag <- mxae / mnae <= convRatio || mxae - mnae <= tol

  # Converged if magnitude and distance are close and error oscillates in sign
  isOscil(errs) && errDistance && errMag
}

# Check denominator polynomial for zero in the requested range
checkDenom <- function(a, l, u) {
  dngrRt <- tryCatch(uniroot(polyCalc, c(l, u), a = a),
                     error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(dngrRt, "simpleError")) {
    return(NULL)
  } else {
    return(dngrRt$root)
  }
}
