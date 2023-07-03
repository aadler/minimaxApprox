# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+


remRatMat <- function(x, E, y, nD, dD) {
  n <- length(x)
  altSgn <- (-1) ^ (seq_len(n) - 1)
  Ev <- altSgn * E
  yv <- -(y + Ev)
  AMat <- vanderMat(x, nD)
  BMat <- vanderMat(x, dD)[, -1] * yv
  cbind(AMat, BMat, -altSgn, deparse.level = 0)
}

remRatCoeffs <- function(x, E, fn, nD, dD) {
  y <- callFun(fn, x)
  P <- solve(remRatMat(x, E, y, nD, dD), y)
  RR <- list(a = P[1:(nD + 1L)],
             b = c(1, P[(nD + 2L):(nD + dD + 1L)]),
             E = P[length(P)])
  if (sum(sapply(RR, length)) != (length(x) + 1)) {
    stop("Catastrophic Error. Result vector not of required length.")
  }
  RR
}

remRatFunc <- function(x, a, b) {
  sum(a * x ^ (seq_along(a) - 1)) / sum(b * x ^ (seq_along(b) - 1))
}

remRatErr <- function(x, a, b, fn) {remRatFunc(x, a, b) - callFun(fn, x)}

remRatRoots <- function(x, a, b, fn) {
  r <- double(length(x) - 1)
  for (i in seq_along(r)) {
    r[i] <- uniroot(remRatErr, lower = x[i], upper = x[i + 1L],
                    a = a, b = b, fn = fn)$root
  }
  r
}

remRatSwitch <- function(r, l, u, a, b, fn, eps) {
  nodes <- data.frame(lower = c(l, r), upper = c(r, u))
  x <- double(length(r) + 1)
  for (i in seq_along(x)) {
    x[i] <- optimize(remRatErr,
                     lower = nodes$lower[i],
                     upper = nodes$upper[i],
                     a = a, b = b, fn = fn,
                     maximum = remRatErr(nodes$lower[i] + eps, a = a, b = b,
                                         fn = fn) > 0)[[1]]
  }
  x
}

remRat <- function(fn, lower, upper, numerd, denomd, opts = list()) {

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
    tol <- 1e-12
  }

  if ("eps" %in% names(opts)) {
    tol <- opts$eps
  } else {
    eps <- min((upper - lower) / 1000, 1e-3)
  }

  n <- numerd + denomd + 2
  x <- chebNodes(n, lower, upper)

  convergeErr <- function(x, fn, tol, nD, dD) {
    E <- 1
    j <- 0
    repeat {
      j <- j + 1
      if (j > maxiter / 100) break
      RR <- remRatCoeffs(x, E, fn, nD, dD)
      if (isTRUE(all.equal(abs(E), abs(RR$E), tol = tol))) break
      E <- (E + RR$E) / 2
    }
    RR
  }

  RR <- convergeErr(x, fn, tol, numerd, denomd)
  i <- 0L
  repeat {
    if (i > maxiter) {
      message("Convergence not acheived in ", maxiter, " iterations.")
      break
    }

    i <- i + 1L
    r <- remRatRoots(x, RR$a, RR$b, fn)
    x <- remRatSwitch(r, lower, upper, RR$a, RR$b, fn, eps)
    RR <- convergeErr(x, fn, tol, numerd, denomd)

    errs <- sapply(x, remRatErr, a = RR$a, b = RR$b, fn = fn)
    abserrs <- abs(errs)
    mxae <- max(abs(errs))
    if (showProgress) message("i: ", i, " E: ", RR$E, " maxErr: ", max(abserrs),
                              " Diff:", max(abserrs) - abs(RR$E))

    if (all(diff(abserrs) < tol) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= abs(RR$E) ||
         isTRUE(all.equal(max(abserrs), abs(RR$E), tol = tol)))) break
  }
  ret <- list(a = RR$a, b = RR$b, EE = RR$E, OE = mxae, iterations = i,
              basis = x)
  attr(ret, "type") <- "Rational"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  class(ret) <- c("RatApprox", class(ret))
  ret
}

