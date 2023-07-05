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
  # P <- solve(remRatMat(x, E, y, nD, dD), y)
  PQ <- qr(remRatMat(x, E, y, nD, dD), LAPACK = TRUE)
  P <- solve.qr(PQ, y)
  RR <- list(a = P[seq_len(nD + 1)],            # Works even if nD = 0
             b = c(1, P[seq_len(dD) + nD + 1]), # Works even if dD = 0
             E = P[length(P)])
  if (sum(lengths(RR)) != (length(x) + 1)) {
    stop("Catastrophic Error. Result vector not of required length.")
  }
  RR
}

# Function to calculate value of minimax rational approximation at x given a & b
remRatFunc <- function(x, a, b)  polyCalc(x, a) / polyCalc(x, b)

# Function to calculate error between known and calculated values
remRatErr <- function(x, a, b, fn) remRatFunc(x, a, b) - callFun(fn, x)

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima
remRatRoots <- function(x, a, b, fn) {
  if (all(abs(remRatErr(x, a, b, fn)) <= 5 * .Machine$double.eps)) {
    stop("This code only functions to machine double precision. All error ",
         "values are below machine double precision. Please try again using a ",
         "lesser degree.")
  }
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remRatErr, interval = intv, a = a, b = b, fn = fn),
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
remRatSwitch <- function(r, l, u, a, b, fn) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  maximize <- sign(remRatErr(l, a, b, fn)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    extrma <- tryCatch(optimize(remRatErr, interval = intv,
                                a = a, b = b, fn = fn, maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    if (inherits(extrma, "simpleError")) {
      endPtErr <- remRatErr(intv, a, b, fn)
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
    E <- remRatErr(p, a, b, fn)

    if (maximize) {
      x[i] <- p[which.max(E)]
    } else {
      x[i] <- p[which.min(E)]
    }

    # Flip maximize
    maximize <- !maximize
  }

  # Test Oscillation
  # if (!isOscil(remRatErr(x, a, b, fn))) {
  #   stop("Control points do not result in oscillating errors.")
  # }
  x
}

# Main function to calculate and return the minimax rational approximation
remRat <- function(fn, lower, upper, numerd, denomd, xi = NULL, opts = list()) {

  # Handle configuration options
  if ("maxiter" %in% names(opts)) {
    maxiter <- opts$maxiter
  } else {
    maxiter <- 3000L
  }

  if ("miniter" %in% names(opts)) {
    miniter <- opts$miniter
  } else {
    miniter <- 15L
  }

  if ("showProgress" %in% names(opts)) {
    showProgress <- opts$showProgress
  } else {
    showProgress <- FALSE
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- 1e-10
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
  # equations until E converges. This function has no other use so defined
  # INSIDE the remRat call.
  convergeErr <- function(x, fn, tol, nD, dD) {
    E <- 1
    j <- 0L
    repeat {
      j <- j + 1L
      if (j > max(maxiter / 100, 1)) break  # Otherwise this would take FOREVER
      RR <- remRatCoeffs(x, E, fn, nD, dD)
      if (all(abs(abs(E) - abs(RR$E)) < tol)) break
      E <- (E + RR$E) / 2
    }
    RR
  }

  RR <- convergeErr(x, fn, tol, numerd, denomd)
  i <- 0L
  repeat {
    i <- i + 1L
    r <- remRatRoots(x, RR$a, RR$b, fn)
    x <- remRatSwitch(r, lower, upper, RR$a, RR$b, fn)
    # if (!is.null(dngr)) {
    #   xrat <- x / dngr
    #   closebelow <- which(xrat <= 1 & xrat > 0.99)
    #   closeabove <- which(xrat > 1 & xrat <= 1.01)
    #   x[closebelow] <- x[closebelow] * 0.9
    #   x[closeabove] <- x[closeabove] * 1.1
    #   x <- sort(x)
    # }
    RR <- convergeErr(x, fn, tol, numerd, denomd)
    dngr <- checkDenom(RR$b, lower, upper)
    if (!is.null(dngr)) {
     stop("The ", denomd, " degree polynomial in the denominator has a zero ",
           "at ", fN(dngr), " which makes rational approximation perilous for ",
           "this function over the interval [", lower, ", ", upper, "].")
    }
    errs <- remRatErr(x, a = RR$a, b = RR$b, fn = fn)
    mxae <- max(abs(errs))

    if (showProgress) {
      message("i: ", i, " E: ", fN(RR$E), " maxErr: ", fN(mxae),
              " Diff:", fN(abs(mxae - abs(RR$E))))
    }

    if ((isConverged(errs, RR$E, tol) && i >= miniter) || i > maxiter) break
  }

  gotWarning <- FALSE

  if (i >= maxiter) {
    mess <- paste("Convergence not acheived in", maxiter, "iterations.\n")
    mess <- paste0(mess, "Maximum observed error ",
                   formatC(mxae / RR$E, digits = 6L), " times expected.")
    warning(mess)
    gotWarning <- TRUE
  }

  if (all(abs(errs) < 10 * .Machine$double.eps)) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
    gotWarning <- TRUE
  }

  ret <- list(a = RR$a, b = RR$b, EE = abs(RR$E), OE = mxae,
              Diff = abs(abs(RR$E) - mxae), iterations = i, x = x,
              Warning = gotWarning)
  attr(ret, "type") <- "Rational"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  class(ret) <- c("RatApprox", class(ret))
  ret
}
