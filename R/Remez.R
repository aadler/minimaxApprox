# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

chebNodes <- function(n, a, b) {
  zapsmall(
    sort(0.5 * (a + b + (b - a) * cos((2 * seq_len(n) - 1) * pi / (2 * n))))
  )
}

vanderMat <- function(x, n) {
  n1 <- n + 1L
  matrix(rep(x, each = n1) ^ (seq_len(n1) - 1L), ncol = n1, byrow = TRUE)
}

callFun <- function(fn, x) {
  if (is.function(fn)) {
    return(do.call(match.fun(fn), args = list(x = x)))
  } else {
    stop("Unable to parse function.")
  }
}

isOscil <- function(x) all(abs(diff(sign(x))) == 2)

polyCalc <- function(x, a) sum(a * x ^ (seq_along(a) - 1L))

remPolyMat <- function(x) {
  n <- length(x)
  A <- vanderMat(x, n - 2)
  cbind(A, (-1) ^ (seq_len(n) - 1))
}

remPolyCoeffs <- function(x, fn) {
  PP <- solve(remPolyMat(x), callFun(fn, x))
  list(b = PP[-length(PP)], E = PP[length(PP)])
}

remPolyErr <- function(x, b, fn) polyCalc(x, b) - callFun(fn, x)

remPolyRoots <- function(x, b, fn) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    r[i] <- uniroot(remPolyErr, lower = x[i], upper = x[i + 1L],
                    b = b, fn = fn)$root
  }
  r
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
    tol <- 1e-12
  }

  # Initial x's
  x <- chebNodes(degree + 2, lower, upper)

  # Initial Polynomial Guess
  PP <- remPolyCoeffs(x, fn)
  i <- 0L
  # Remez Loop. Must be performed at least once so use "do-until" logic
  repeat {
    if (i > maxiter) {
      stop("Convergence not acheived in ", maxiter, " iterations.")
    }
    i <- i + 1L
    r <- remPolyRoots(x, PP$b, fn)
    x <- remPolySwitch(r, lower, upper, PP$b, fn)
    PP <- remPolyCoeffs(x, fn)
    errs <- sapply(x, remPolyErr, b = PP$b, fn = fn)
    abserrs <- abs(errs)
    if (showProgress) message(paste("i:", i, "E:", PP$E, "maxErr:", max(abserrs)))

    if (all(diff(abserrs) < tol) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= abs(PP$E) ||
         isTRUE(all.equal(max(abserrs), abs(PP$E), tol = tol)))) break
  }
  ret <- list(b = PP$b, E = PP$E)
  attr(ret, "type") <- "Polynomial"
  attr(ret, "iterations") <- i
  attr(ret, "basis") <- x
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  class(ret) <- c("RatApprox", class(ret))

  ret
}

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
    if (showProgress) message("i: ", i, " E: ", RR$E, " maxErr: ", max(abserrs),
                              " Diff:", max(abserrs) - abs(RR$E))

    if (all(diff(abserrs) < tol) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= abs(RR$E) ||
         isTRUE(all.equal(max(abserrs), abs(RR$E), tol = tol)))) break
  }
  ret <- list(a = RR$a, b = RR$b, E = RR$E)
  attr(ret, "type") <- "Rational"
  attr(ret, "iterations") <- i
  attr(ret, "basis") <- x
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  class(ret) <- c("RatApprox", class(ret))
  ret
}

print.RatApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    print(list(b = x$b, E = x$E))
  } else {
    print(list(a = x$a, b = x$b, E = x$E))
  }
}

plot.RatApprox <- function(x, ...) {
  rng <- attr(x, "range")
  xF <- attr(x, "basis")
  z <- seq(rng[1], rng[2], length.out = 1000L)

  if (attr(x, "type") == "Polynomial") {
    zz <- sapply(z, remPolyErr, x$b, attr(x, "func"))
    yF <- sapply(xF, remPolyErr, x$b, attr(x, "func"))
  } else {
    zz <- sapply(z, remRatErr, x$a, x$b, attr(x, "func"))
    yF <- sapply(xF, remRatErr, x$a, x$b, attr(x, "func"))
  }

  plot(z, zz, type = 'l',  xlab = "x", ylab = "Error")
  abline(h = 0)
  points(xF, yF, col = "red", pch = 16)
}
