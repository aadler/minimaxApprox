# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, errType = "abs", xi = NULL,
                          opts = list()) {

  absErr <- switch(tolower(substring(errType, 1, 3)),
                   abs = TRUE,
                   rel = FALSE,
                   stop("Error type must be either 'abs'olute or 'rel'ative."))

  if (length(degree) == 2L) {        # If rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
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
    remPoly(fn, lower, upper, as.integer(degree), absErr, opts)
  }
}

# Evaluation convenience function
minimaxEval <- function(x, mmA) {
  if (attr(mmA, "type") == "Polynomial") {
    polyCalc(x, mmA$a)
  } else {
    remRatFunc(x, mmA$a, mmA$b)
  }
}

# Print method (hide i and basis/x but leave in list and not in attribute)
print.minimaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    coefficients <- list(a = x$a)
  } else {
    coefficients <- list(a = x$a, b = x$b)
  }

  diags <- list(x$EE,
                x$OE,
                Ratio = round(x$OE / x$EE, 6L),
                Difference = abs(x$OE - x$EE),
                Warnings = x$Warning)

  names(diags)[1:2] <- if (attr(x, "absErr")) {
     c("ExpectedAbsError", "ObservedAbsError")
  } else {
    c("ExpectedRelError", "ObservedRelError")
  }

  print(c(coefficients, diags))
}

coef.minimaxApprox <- function(object, ...) {
  if (attr(object, "type") == "Polynomial") {
    coef <- list(a = object$a)
  } else {
    coef <- list(a = object$a, b = object$b)
  }

  coef
}
# Plot method for errors and basis points
plot.minimaxApprox <- function(x, y, ...) {
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

  plot(z, zz, type = "l", xlab = "x", ylab = ylab, ...)
  abline(h = 0)
  points(x$x, y, col = "red", pch = 16L)
  abline(h = c(-x$OE, x$OE), lty = 2L, col = "blue")
  abline(h = c(-x$EE, x$EE), lty = 3L, col = "red")
  legend(x = "bottomleft", inset = c(0.35, 1), col = c("red", "red", "blue"),
         lty = c(NA, 3L, 2L), legend = c("Basis", "Exp Err", "Obs Err"),
         pch = c(16L, NA, NA), bg = "transparent", xpd = TRUE)
}
