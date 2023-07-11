# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, errType = "abs", xi = NULL,
                          opts = list()) {

  absErr <- switch(tolower(substring(errType, 1, 3)),
                   abs = TRUE,
                   rel = FALSE,
                   stop("Error type must be either 'abs'olute or 'rel'ative."))

  if (length(degree) == 2L) {        # If rational approximation requested
    numerd <- degree[1L]
    denomd <- degree[2L]
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE               #  Polynomial approximation requested
    if (!is.null(xi)) {
      warning("Polynomial approximation uses Chebyeshev nodes for initial ",
              "guess. Any passed xi is ignored.")
    }
  } else {
    # All else is an error
    stop("Polynomial approximation takes one value for degree and rational ",
         "approximation takes a vector of two values for numerator and ",
         "denominator. Any other inputs are invalid.")
  }

  if (ratApprox) {
    remRat(fn, lower, upper, numerd, denomd, absErr, xi, opts)
  } else {
    remPoly(fn, lower, upper, degree, absErr, opts)
  }
}

# Print method (hide i and basis/x but leave in list and not in attribute)
print.minimaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    ret <- list(a = x$a)
  } else {
    ret <- list(a = x$a, b = x$b)
  }

  ret <- c(ret, list(x$EE,
                     x$OE,
                     Ratio = round(x$OE / x$EE, 6L),
                     Difference = abs(x$OE - x$EE),
                     Warnings = x$Warning))
  if (attr(x, "absErr")) {
    names(ret)[3:4] <- c("ExpectedAbsError", "ObservedAbsError")
  } else {
    names(ret)[3:4] <- c("ExpectedRelError", "ObservedRelError")
  }
  print(ret)
}

coef.minimaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    coef <- list(a = x$a)
  } else {
    coef <- list(a = x$a, b = x$b)
  }

  coef
}
# Plot method for errors and basis points
plot.minimaxApprox <- function(x, ...) {
  args <- list(...)
  rng <- attr(x, "range")
  fn <- attr(x, "func")
  absErr <- attr(x, "absErr")
  z <- seq(rng[1], rng[2], length.out = 1001L)

  if (attr(x, "type") == "Polynomial") {
    zz <- remPolyErr(z, x$a, fn, absErr)
    y <- remPolyErr(x$x, x$a, fn, absErr)
  } else {
    zz <- remRatErr(z, x$a, x$b, fn, absErr)
    y <- remRatErr(x$x, x$a, x$b, fn, absErr)
  }

  ylab <- if (absErr) "Absolute Error" else "Relative Error"
  ybnd <- max(abs(x$EE), abs(x$OE))
  if ("ylim" %in% names(args)) {
    ylim <- args$ylim
  } else {
    ylim <- c(-ybnd, ybnd)
  }

  plot(z, zz, type = "l", xlab = "x", ylab = ylab, ylim = ylim)
  abline(h = 0)
  points(x$x, y, col = "red", pch = 16)
  abline(h = c(-x$EE, x$EE), lty = 2, col = "red")
  abline(h = c(-x$OE, x$OE), lty = 3, col = "blue")
  legend(x = "bottomleft", inset = c(0.35, 1), col = c("red", "red", "blue"),
         lty = c(NA, 2, 3), legend = c("Basis", "Exp Err", "Obs Err"),
         pch = c(16, NA, NA), bg = "transparent", xpd = TRUE)
}
