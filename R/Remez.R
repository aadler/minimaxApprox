chebNodes <- function(n, a, b) {
  zapsmall(
    sort(0.5 * (a + b + (b - a) * cos((2 * seq_len(n) - 1) * pi / (2 * n))))
  )
}

vanderMat <- function(x, n) {
  n1 <- n + 1
  matrix(rep(x, each = n1) ^ (seq_len(n1) - 1), ncol = n1, byrow = TRUE)
}

callFun <- function(fn, x) {
  if (is.character(fn)) {
    y <- do.call(match.fun(eval(parse(text = fn))), args = list(x = x))
  } else if (is.function(fn)) {
    y <- do.call(match.fun(fn), args = list(x = x))
  } else {
    stop("Unable to parse function")
  }
  y
}

remPolyMat <- function(x) {
  n <- length(x)
  A <- vanderMat(x, n - 2)
  cbind(A, (-1) ^ (seq_len(n) - 1))
}

remPolyCoeffs <- function(x, fn) {
  PP <- solve(remPolyMat(x), callFun(fn, x))
  list(b = PP[-length(PP)], E = abs(PP[length(PP)]))
}

remPolyErr <- function(x, b, fn) {
  sum(b * x ^ (seq_along(b) - 1)) - callFun(fn, x)
}

remPolyRoots <- function(x, b, fn) {
  r <- double(length(x) - 1)
  for (i in seq_along(r)) {
    r[i] <- uniroot(remPolyErr, lower = x[i], upper = x[i + 1L],
                    b = b, fn = fn)$root
  }
  r
}

remPolySwitch <- function(r, l, u, b, fn, eps) {
  nodes <- data.frame(lower = c(l, r), upper = c(r, u))
  x <- double(length(r) + 1)
  for (i in seq_along(x)) {
    x[i] <- optimize(remPolyErr,
                     lower = nodes$lower[i],
                     upper = nodes$upper[i],
                     b = b, fn = fn,
                     maximum = remPolyErr(nodes$lower[i] + eps,
                                      b = b, fn = fn) > 0)[[1]]
  }
  x
}

remPoly <- function(fn, lower, upper, degree, showProgress = FALSE, opts = list()) {

  # Handle configuration options
  if ("maxiter" %in% names(opts)) {
    maxiter <- opts$maxit
  } else {
    maxiter <- 2500
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- sqrt(.Machine$double.eps)
  }

  if ("eps" %in% names(opts)) {
    tol <- opts$eps
  } else {
    eps <- min((upper - lower) / 1000, 1e-3)
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
    x <- remPolySwitch(r, lower, upper, PP$b, fn, eps)
    PP <- remPolyCoeffs(x, fn)
    errs <- sapply(x, remPolyErr, b = PP$b, fn = fn)
    abserrs <- abs(errs)
    if (showProgress) message(paste("i:", i, "E:", PP$E, "maxErr:", max(abserrs)))

    if (all(diff(abserrs) < tol) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= PP$E ||
         isTRUE(all.equal(max(abserrs), PP$E, tol = tol)))) break
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
             E = abs(P[length(P)]))
  if (sum(sapply(RR, length)) != (length(x) + 1)) {
    stop("Catastrophic Error. Result vector not of required length.")
  }
  RR
}

remRatErr <- function(x, a, b, fn) {
  sum(a * x ^ (seq_along(a) - 1)) / sum(b * x ^ (seq_along(b) - 1)) -
    callFun(fn, x)
}

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

remRat <- function(fn, lower, upper, numerd, denomd, showProgress = FALSE,
                   opts = list()) {

  # Handle configuration options
  if ("maxiter" %in% names(opts)) {
    maxiter <- opts$maxit
  } else {
    maxiter <- 2500
  }

  if ("tol" %in% names(opts)) {
    tol <- opts$tol
  } else {
    tol <- sqrt(.Machine$double.eps)
  }

  if ("eps" %in% names(opts)) {
    tol <- opts$eps
  } else {
    eps <- min((upper - lower) / 1000, 1e-3)
  }

  n <- numerd + denomd + 2
  x <- chebNodes(n, lower, upper)

  convergeErr <- function(x, fn, tol, nD, dD) {
    E <- 0
    j <- 0
    repeat {
      j <- j + 1L
      RR <- remRatCoeffs(x, E, fn, nD, dD)
      if (isTRUE(all.equal(abs(E), RR$E, tol = tol))) break
      E <- (abs(E) + RR$E) / 2
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
                              " Diff:", max(abserrs) - RR$E)

    if (all(diff(abserrs) < tol) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= RR$E ||
         isTRUE(all.equal(max(abserrs), RR$E, tol = tol)))) break
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
  z <- seq(rng[1], rng[2], length.out = 1000L)

  if (attr(x, "type") == "Polynomial") {
    zz <- sapply(z, remPolyErr, x$b, attr(x, "func"))
  } else {
    zz <- sapply(z, remRatErr, x$a, x$b, attr(x, "func"))
  }
  plot(z, zz, type = 'l',  xlab = "x", ylab = "Error")
  abline(h = 0)
  abline(v = attr(x, "basis"), col = "red")
}
