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

remPoly <- function(fn, lower, upper, degree, eps = 1e-3,
                    tol = .Machine$double.eps, showProgress = FALSE) {
  # Initial x's
  x <- chebNodes(degree + 2, lower, upper)
  # Initial Polynomial Guess
  PP <- remPolyCoeffs(x, fn)
  i <- 0L
  # Remez Loop. Must be performed at least once so use "do-until" logic
  repeat {
    i <- i + 1L
    r <- remPolyRoots(x, PP$b, fn)
    x <- remPolySwitch(r, lower, upper, PP$b, fn, eps)
    PP <- remPolyCoeffs(x, fn)
    errs <- sapply(x, remPolyErr, b = PP$b, fn = fn)
    abserrs <- abs(errs)
    if (showProgress) message(paste("i:", i, "E:", PP$E, "maxErr:", max(abserrs)))
    # Run so long as max(abserr) > E
    if (all(diff(abserrs) < .Machine$double.eps) &&
        all(abs(diff(sign(errs))) == 2) &&
        (max(abserrs) <= PP$E || isTRUE(all.equal(max(abserrs), PP$E)))) break
  }
  return(list(b = PP$b, E = PP$E, iter = i, x = x))
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

remRat <- function(fn, lower, upper, numerd, denomd, eps = 1e-5,
                   tol = .Machine$double.eps, showProgress = FALSE) {
  n <- numerd + denomd + 2
  x <- chebNodes(n, lower, upper)
  y <- callFun(fn, x)

  convergeErr <- function(x, fn, tol, nD, dD) {
    E <- 1
    repeat {
      RR <- remRatCoeffs(x, E, fn, nD, dD)
      if (isTRUE(all.equal(E, RR$E, tol = tol))) break
      E <- RR$E
    }
    RR
  }

  RR <- convergeErr(x, fn, tol, numerd, denomd)
  i <- 0L
  repeat {
    i <- 1 + 1L
    r <- remRatRoots(x, RR$a, RR$b, fn)
    x <- remRatSwitch(r, lower, upper, RR$a, RR$b, fn, eps)
    RR <- convergeErr(x, fn, tol, numerd, denomd)
    abserr <- abs(sapply(x, remRatErr, a = RR$a, b = RR$b, fn = fn))
    if (showProgress) message(paste("i:", i, "E:", RR$E, "maxErr:",
                                    max(abserr), "Diff:", max(abserr) - RR$E))
    # Run so long as max(abserr) > E
    if (max(abserr) <= RR$E || isTRUE(all.equal(max(abserr), RR$E))) break
  }
  RR
}
