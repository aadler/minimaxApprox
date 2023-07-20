# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, relErr = FALSE, xi = NULL,
                          opts = list()) {

  # Handle configuration options
  if (!("maxiter" %in% names(opts))) {
    opts$maxiter <- 100L
  }

  if (!("miniter" %in% names(opts))) {
    opts$miniter <- 10L
  }

  if (!("conviter" %in% names(opts))) {
    opts$conviter <- 10L
  } else {
    # If actually passed then overwrite maxiter if conviter > maxiter
    opts$maxiter <- max(opts$maxiter, opts$conviter)
  }

  if (!("showProgress" %in% names(opts))) {
    opts$showProgress <- FALSE
  }

  if (!("convrat" %in% names(opts))) {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    opts$convrat <- 1.000000001
  }

  if (!("tol" %in% names(opts))) {
    opts$tol <- 1e-14
  }

  if (!is.logical(relErr)) {
    stop("Relative Error must be a logical value. ",
         "Default FALSE returns absolute error.")
  }

  if (length(degree) == 2L) {        # Rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE               # Polynomial approximation requested
    if (!is.null(xi)) {
      message("Polynomial approximation uses Chebyshev nodes for initial ",
              "guess. Any passed xi is ignored.")
    }
  } else {
    # All else is an error
    stop("Polynomial approximation takes one value for degree and rational ",
         "approximation takes a vector of two values for numerator and ",
         "denominator degrees. Any other inputs are invalid.")
  }

  # Call Calculation Functions
  mmA <- if (ratApprox) {
    remRat(fn, lower, upper, numerd, denomd, relErr, xi, opts)
  } else {
    remPoly(fn, lower, upper, as.integer(degree), relErr, opts)
  }

  # Handle all warnings centrally
  gotWarning <- FALSE

  if (mmA$i >= opts$maxiter && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not achieved in ",
            mmA$i, " iterations.\n", "The ratio is ", fC(mmA$mxae / mmA$expe),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$unchanged && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not achieved in ",
            mmA$i, " iterations.\n", mmA$unchanging_i, " successive calculated ",
            "errors were too close to each other to warrant further ",
            "iterations.\nThe ratio is ", fC(mmA$mxae / mmA$expe), " times ",
            "expected and the difference is ", fC(abs(mmA$mxae - mmA$expe)),
            " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$mxae < 10 * .Machine$double.eps) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
    gotWarning <- TRUE
  }

  coefficients <- if (ratApprox) list(a = mmA$a, b = mmA$b) else list(a = mmA$a)
  diagnostics <- list(EE = mmA$expe, OE = mmA$mxae,  iterations = mmA$i,
                      x = mmA$x, Warning = gotWarning)
  ret <- c(coefficients, diagnostics)
  attr(ret, "type") <- if (ratApprox) "Rational" else "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "relErr") <- relErr
  attr(ret, "tol") <- opts$tol
  attr(ret, "convrat") <- opts$convrat
  class(ret) <- c("minimaxApprox", class(ret))

  ret
}

# Evaluation convenience function. Identical to evalFunc but tests for
# inheritance from minimaxApprox.
minimaxEval <- function(x, mmA) {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  evalFunc(x, mmA)
}
