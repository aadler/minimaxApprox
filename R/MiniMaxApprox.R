# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

minimaxApprox <- function(fn, lower, upper, degree, relErr = FALSE,
                          basis = "Chebyshev", xi = NULL, opts = list()) {

  ## ------------------------------------------------------------------------
  ## Input validation (F2/F10). This entire block runs before any option
  ## default is set or option-derived quantity is computed. In particular,
  ## opts$tailtol's default below uses (upper - lower), so lower/upper must
  ## be validated first, or a bad range would silently poison that default
  ## before ever being checked itself.
  ## ------------------------------------------------------------------------

  # fn: must be a function whose first formal argument is 'x' -- the
  # documented contract (man/MiniMaxApprox.Rd), required because callFun()
  # always dispatches via do.call(fn, list(x = x)). Primitives (exp, sin,
  # gamma, ...) report formals(fn) == NULL, since their argument list is not
  # an R closure; args(fn) returns a stub closure with a real formals list
  # for primitives, and is a harmless no-op for ordinary closures -- so it is
  # used uniformly here instead of formals(fn) directly, which would
  # incorrectly reject every primitive.
  if (!is.function(fn)) {
    stop("'fn' must be a function whose first argument is 'x'.")
  }
  fnArgs <- names(formals(args(fn)))
  if (length(fnArgs) == 0L || fnArgs[1L] != "x") {
    stop("'fn' must be a function whose first argument is 'x'.")
  }

  # lower/upper: finite, non-missing numeric scalars, with lower < upper.
  # (F2: an inverted range was previously never checked and silently
  # returned a badly suboptimal result instead of erroring.)
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower)) {
    stop("'lower' must be a finite, non-missing numeric scalar.")
  }
  if (!is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("'upper' must be a finite, non-missing numeric scalar.")
  }
  if (lower >= upper) {
    stop("'lower' must be less than 'upper'. Did you mean to swap the ",
         "arguments?")
  }

  # degree: finite, non-missing numeric. Length (1 vs 2 vs invalid) is still
  # dispatched further down with its existing message; this guard only
  # ensures degree is safe to compare/floor -- previously, e.g.,
  # `any(NA < 0)` is NA, and `if (NA)` threw "missing value where TRUE/FALSE
  # needed" instead of a real error.
  if (!is.numeric(degree) || anyNA(degree) || !all(is.finite(degree))) {
    stop("'degree' must be finite, non-missing numeric value(s).")
  }

  # relErr (existing check, unchanged; relocated into the validation block).
  if (!is.logical(relErr)) {
    stop("Relative Error must be a logical value. ",
         "Default FALSE returns absolute error.")
  }

  # basis (existing check, unchanged; relocated into the validation block).
  basis <- tolower(substr(basis, 1L, 1L))
  if (!(basis %in% c("c", "m"))) {
    stop("Must select either 'C'hebyshev or 'm'onomial basis for analysis.")
  }

  # opts: must be a list; only members the caller actually supplied are
  # validated here (unknown names are still silently ignored downstream,
  # unchanged behavior). This runs before the defaults block below, so a bad
  # value (e.g. maxiter = 0) cannot slip through -- the defaults block below
  # only fills in *missing* names, it never overwrites a supplied one.
  if (!is.list(opts)) {
    stop("'opts' must be a list.")
  }
  nopts <- names(opts)

  chkPosInt <- function(val, nm) {
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        !is.finite(val) || val < 1 || floor(val) != val) {
      stop("'opts$", nm, "' must be a single positive integer.")
    }
  }
  chkNumScalar <- function(val, nm) {
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        !is.finite(val)) {
      stop("'opts$", nm, "' must be a single finite, non-missing numeric ",
           "value.")
    }
  }
  for (nm in c("maxiter", "miniter", "conviter")) {
    if (nm %in% nopts) chkPosInt(opts[[nm]], nm)
  }
  for (nm in c("tol", "convrat")) {
    if (nm %in% nopts) chkNumScalar(opts[[nm]], nm)
  }
  # tailtol/ztol may deliberately be NULL by design (NULL disables the
  # tailtol restart check / requests no coefficient zeroing); only validate
  # when the caller supplied a non-NULL value.
  for (nm in c("tailtol", "ztol")) {
    if (nm %in% nopts && !is.null(opts[[nm]])) chkNumScalar(opts[[nm]], nm)
  }

  ## ------------------------------------------------------------------------
  ## End input validation.
  ## ------------------------------------------------------------------------

  # Handle configuration options
  if (!("maxiter" %in% nopts)) {
    opts$maxiter <- 100L
  }

  if (!("miniter" %in% nopts)) {
    opts$miniter <- 10L
  }

  if ("conviter" %in% nopts) {
    # If actually passed then overwrite both maxiter and miniter if conviter is
    # greater than either one.
    opts$maxiter <- max(opts$maxiter, opts$conviter)
    opts$miniter <- max(opts$miniter, opts$conviter)
  } else {
    opts$conviter <- 30L
  }

  if (!("showProgress" %in% nopts)) {
    opts$showProgress <- FALSE
  }

  if (!("convrat" %in% nopts)) {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    opts$convrat <- 1.000000001
  }

  if (!("tol" %in% nopts)) {
    opts$tol <- 1e-14
  }

  # Used for cases where we check polynomial degree n + 1.
  # See issue 2 https://github.com/aadler/minimaxApprox/issues/2
  if (!("tailtol" %in% nopts)) {
    opts$tailtol <- min(1e-10, (upper - lower) / 1e6)
  }

  if (!("ztol" %in% nopts)) {
    opts$ztol <- NULL
  }

  if (any(degree < 0) || any(floor(degree) < degree)) {
    stop("Degrees must be integers of least 0 (constant).")
  }

  if (length(degree) == 2L) {         # Rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE                # Polynomial approximation requested
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
    remRat(fn, lower, upper, numerd, denomd, relErr, basis, xi, opts)
  } else {
    tryCatch(remPoly(fn, lower, upper, as.integer(degree), relErr, basis, opts),
             error = function(e) simpleError(trimws(e$message)))
  }

  # In response to issue 2, https://github.com/aadler/minimaxApprox/issues/2,
  # allow the polynomial algorithm to try degree n + 1 if it fails due to
  # singular error in degree n. IF the resulting highest coefficient contributes
  # less than opts$tailtol to the result then consider it 0 and return the
  # resulting degree n and message appropriately. For rational approximation or
  # if the error is not a simpleError or does not contain the word "singular",
  # let the default R failure message come through.
  #
  # TODO: Trap other errors with better messages (AA: 2025-12-24)


  if (!ratApprox && inherits(mmA, "simpleError") &&
      grepl("singular", mmA$message, fixed = TRUE)) {

    if (is.null(opts$tailtol)) {
      stop("The algorithm did not converge when looking for a polynomial of ",
           "degree ", degree, " and NULL was passed to the tailtol option.")
    }

    mmA <- tryCatch(remPoly(fn, lower, upper, as.integer(degree + 1L),
                            relErr, basis, opts),
                    error = function(e) simpleError(trimws(e$message)))

    if (inherits(mmA, "simpleError")) {
      # F4: both the degree-n and degree-(n+1) Remez solves failed singular.
      # Before raising the original hard error, check whether the "failure" is
      # actually that fn is (near-)exactly a polynomial of degree <= n, i.e.
      # there is no minimax problem left to solve. interpRescue() builds the
      # plain degree-n interpolant (NOT a Remez result) and returns it only if
      # its error is at the machine-precision floor; otherwise NULL. If NULL,
      # the singularity had a genuine (non-representability) cause and the
      # original error is raised unchanged -- the rescue never masks it.
      rescue <- interpRescue(fn, lower, upper, as.integer(degree), relErr,
                             basis)

      if (is.null(rescue)) {
        stop("The algorithm neither converged when looking for a polynomial of",
             " degree ", degree, " nor when looking for a polynomial of degree ",
             degree + 1L, ".")
      }

      # Assemble an mmA-like result from the interpolant. expe and mxae are
      # both the probe error (there is no leveled E from a solve; ratio is 1).
      # x is set to the interpolation nodes (the interpolant's natural
      # reference; these are NOT equioscillation extrema). The `rescued` flag
      # triggers the not-a-Remez warning in the central warning block below.
      mmA <- list(a = rescue$a, expe = rescue$err, mxae = rescue$err,
                  i = 0L, x = rescue$x, converged = TRUE, unchanged = FALSE,
                  unchanging_i = 0L, zeroBasisError = FALSE, rescued = TRUE)

    } else {
      # Degree-(n+1) retry succeeded: existing tailtol "uppermost coefficient
      # effectively 0" logic, unchanged.
      n <- length(mmA$a)
      # F3 fix: was (mmA$a[n] * xmax^(n-1L)) > opts$tailtol -- no abs() on the
      # coefficient, so any NEGATIVE top coefficient of arbitrary magnitude
      # passed this test and was silently dropped as "effectively zero". Now
      # routed through tailContribution, which takes abs() and uses the
      # basis-correct scale (Chebyshev's unmapped-basis bound differs from
      # xmax^(n-1); see basisScale in shared.R).
      if (tailContribution(mmA$a[n], n, lower, upper, basis) > opts$tailtol) {
        stop("The algorithm did not converge when looking for a polynomial of",
             " degree ", degree, " and when looking for a polynomial of degree ",
             degree + 1L, " the uppermost coefficient is not effectively zero.")
      }

      mmA$a <- mmA$a[-n]
      message("The algorithm failed while looking for a polynomial of degree ",
              degree, " but successfully completed when looking for a",
              " polynomial of degree ", degree + 1L, " with the largest",
              " coefficient's contribution to the approximation <= the tailtol",
              " option. The result is a polynomial of degree ", degree, " as",
              " the uppermost coefficient is effectively 0.")
    }
  }

  # F4, second manifestation. On some BLAS/LAPACK platforms an
  # exactly-representable / precision-resolved function does NOT make the
  # augmented solve report singular; instead the Remez iteration runs to
  # maxiter (or stalls "unchanging") wandering at the machine-precision floor,
  # and mmA comes back as a completed-but-not-converged result rather than a
  # simpleError. Catch that here, DOWNSTREAM of the entire singular/restart
  # block above (mmA is already a finished object; remPoly, switchX, and the
  # restart machinery have all run), so this cannot interfere with the
  # load-bearing singular->restart path the way an upstream guard would (cf. M3
  # reverted collapse-guard). interpRescue()'s own probe -- including the relErr
  # zero-of-fn guard -- is the sole arbiter: it returns the clean interpolant
  # only if the function really is resolved to the machine floor at degree n,
  # otherwise NULL (leaving a genuine non-convergence untouched, so its normal
  # maxiter/unchanging warning still fires).
  if (!ratApprox && !inherits(mmA, "simpleError") && !isTRUE(mmA$rescued) &&
      !mmA$converged) {
    rescue <- interpRescue(fn, lower, upper, as.integer(degree), relErr, basis)
    if (!is.null(rescue)) {
      mmA <- list(a = rescue$a, expe = rescue$err, mxae = rescue$err,
                  i = 0L, x = rescue$x, converged = TRUE, unchanged = FALSE,
                  unchanging_i = 0L, zeroBasisError = FALSE, rescued = TRUE)
    }
  }

  # Handle all warnings centrally.
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
            mmA$i, " iterations.\n", mmA$unchanging_i, " successive ",
            "calculated solutions were too close to each other to warrant ",
            "further iterations.\nThe ratio is ",
            fC(mmA$mxae / mmA$expe, d = 14L),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (isTRUE(mmA$rescued)) {
    warning("The requested degree resolves the function to within machine ",
            "floating-point precision, so no Remez minimax iteration was ",
            "possible (the augmented system is singular because there is no ",
            "minimax problem left to solve). The returned polynomial is the ",
            "degree-", degree, " interpolant through the Chebyshev reference; ",
            "its maximum error (", fC(mmA$mxae), ") is at the ",
            "machine-precision floor, below the level at which a distinct ",
            "minimax solution can be discerned. This is NOT technically a ",
            "Remez result.")
    gotWarning <- TRUE
  }

  # The rescue warning above already states the machine-precision caveat in
  # rescue-specific terms; suppress the generic near-eps warning in that case so
  # a rescued result raises exactly one, more informative warning rather than
  # two overlapping ones.
  if (mmA$mxae < 10 * .Machine$double.eps && !isTRUE(mmA$rescued)) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal given floating point limitations.")
    gotWarning <- TRUE
  }

  if (mmA$zeroBasisError) {
    warning("During convergence, the algorithm chose basis point(s) where the ",
            "functional value is 0. The basis point was perturbed by 1e-12, ",
            "but consider approximating using absolute---not relative---error.")
    gotWarning <- TRUE
  }

  coeff <- if (ratApprox) {
    list(a = mmA$a, b = mmA$b)
  } else {
    list(a = mmA$a)
  }

  if (basis == "m") {
    monomialEq <- NULL
    polynomalBasis <- "Monomial"
  } else {
    monomialEq <- list(aMono = cheb2mon(mmA$a))
    polynomalBasis <- "Chebyshev"
    if (ratApprox) {
      monomialEq <- c(monomialEq, list(bMono = cheb2mon(mmA$b)))
      monomialEq <- mapply(`/`, monomialEq, monomialEq$bMono[1L],
                           SIMPLIFY = FALSE)
    }
  }

  diagnostics <- list(ExpErr = mmA$expe, ObsErr = mmA$mxae, iterations = mmA$i,
                      Extrema = mmA$x, Warning = gotWarning)
  ret <- c(coeff, monomialEq, diagnostics)
  attr(ret, "type") <- if (ratApprox) "Rational" else "Polynomial"
  attr(ret, "basis") <- polynomalBasis
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "relErr") <- relErr
  attr(ret, "tol") <- opts$tol
  attr(ret, "convrat") <- opts$convrat
  class(ret) <- c("minimaxApprox", class(ret))

  ret
}
