# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

remPolyMat <- function(x) {
  n <- length(x)
  A <- vanderMat(x, n - 2)
  cbind(A, (-1) ^ (seq_len(n) - 1))
}

remPolyCoeffs <- function(x, fn) {
  # PPQ <- qr(remPolyMat(x), LAPACK = TRUE, tol = 1e-15)
  # PP <- qr.solve(PPQ, callFun(fn, x), tol = 1e-15)
  PP <- solve(remPolyMat(x), callFun(fn, x))
  list(b = zapsmall(PP[-length(PP)]), E = PP[length(PP)])
}

remPolyErr <- function(x, b, fn) polyCalc(x, b) - callFun(fn, x)

remPolyRoots <- function(x, b, fn, errs) {
  if (all(errs == 0)) {
    return(x)
  } else {
    r <- double(length(x) - 1L)
    for (i in seq_along(r)) {
      intv <- c(x[i], x[i + 1L])
      root <- tryCatch(uniroot(remPolyErr, interval = intv, b = b, fn = fn),
                       error = function(cond) simpleError(trimws(cond$message)))
      if (inherits(root, "simpleError")) {
        r[i] <- intv[which.min(abs(intv))]
      } else {
        r[i] <- root$root
      }
    }
    return(r)
  }
}

remPolySwitch <- function(r, l, u, b, fn) {
  nodes <- data.frame(lower = c(l, r), upper = c(r, u))
  x <- double(length(r) + 1L)
  maximize <- sign(remPolyErr(l, b, fn)) == 1
  for (i in seq_along(x)) {
    x[i] <- optimize(remPolyErr,
                     lower = nodes$lower[i],
                     upper = nodes$upper[i],
                     b = b, fn = fn,
                     maximum = maximize)[[1L]]
    maximize <- !maximize
  }
  errs <- sapply(x, remPolyErr, b = b, fn = fn)
  if (!isOscil(errs)) {
    stop("Control points are not oscillating in sign.\n", x)
  }
  x
}

remPoly <- function(fn, lower, upper, degree, opts = list()) {

  # Handle configuration options
  if ("maxiter" %in% names(opts)) {
    maxiter <- opts$maxit
  } else {
    maxiter <- 2500
  }

  if ("showProgress" %in% names(opts)) {
    maxiter <- opts$showProgress
  } else {
    showProgress <- FALSE
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- 5 * .Machine$double.eps
  }

  # Initial x's
  x <- chebNodes(degree + 2, lower, upper)

  # Initial Polynomial Guess
  PP <- remPolyCoeffs(x, fn)
  errs <- sapply(x, remPolyErr, b = PP$b, fn = fn)
  i <- 0
  while ((!isConverged(errs, PP$E, tol) && i <= maxiter)) {
    i <- i + 1L
    r <- remPolyRoots(x, PP$b, fn, errs)
    x <- remPolySwitch(r, lower, upper, PP$b, fn)
    PP <- remPolyCoeffs(x, fn)
    errs <- sapply(x, remPolyErr, b = PP$b, fn = fn)
    mxae <- max(abs(errs))
    if (showProgress) {
      message("i: ", i, " E: ", fN(PP$E), " maxErr: ", fN(mxae))
    }
  }

  if (i >= maxiter) {
    mess <- paste("Convergence not acheived in", maxiter, "iterations.\n")
    mess <- paste0(mess, "Maximum observed error ",
                   formatC(mxae / PP$E, digits = 4), " times expected.")
    warning(mess)
  }

  ret <- list(b = PP$b, EE = abs(PP$E), OE = mxae, iterations = i,
              basis = x)
  attr(ret, "type") <- "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  class(ret) <- c("RatApprox", class(ret))

  ret
}
