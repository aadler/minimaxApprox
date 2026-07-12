# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Printing convenience function.
fC <- function(x, d = 6L, f = "g", w = -1L) {
  formatC(x, digits = d, format = f, width = w)
}

# Default Chebyshev nodes.
chebNodes <- function(n, a, b) {
  n <- as.integer(n)
  sort(0.5 * (a + b + (b - a) * cos((2 * seq_len(n) - 1) * pi / (2 * n))))
}

# Call the function being approximated on the points x.
callFun <- function(fn, x) {
  if (!is.function(fn)) stop("Unable to parse function.")
  do.call(match.fun(fn), args = list(x = x))
}

# Check that the values passed are oscillating in sign.
# F8 fix: sign(x) on NaN/NA input propagates NA through all(), which later
# errors inside an if() (isConverged is used as an if() condition). !anyNA()
# short-circuits that to a clean FALSE. A zero error (sign 0) is treated as
# non-oscillating: an exact zero at a reference point is not a magnitude-E
# equioscillation extremum, and it is behaviorally inert for convergence anyway
# (a zero forces mnae = 0 in isConverged, which independently fails the
# Magnitude Test).
isOscil <- function(x) {
  s <- sign(x)
  !anyNA(s) && all(abs(diff(s)) == 2)
}

evalFunc <- function(x, R, basis, l, u) {
  # M6 (F5 Option A): l/u are now REQUIRED (not defaulted) on purpose. This is
  # the single innermost dispatch site for Chebyshev-basis evaluation; a
  # missed CHEBYSHEV call site anywhere upstream (findRoots/remErr/plot/
  # minimaxEval/...) now surfaces immediately as "argument l is missing"
  # rather than silently evaluating T_k at raw x. NOTE: R's lazy evaluation
  # means this only fires when l/u are actually dereferenced -- monomial-
  # basis callers that omit l/u run fine (harmless, since this branch is
  # never reached for them), which is why monomial callers below still pass
  # them: for signature consistency with the required-parameter contract,
  # not because omitting them would break anything for basis = "m". See
  # switchX/checkDenom/basisScale for the established precedent of required
  # (not optional) l/u on internals that need the range.
  z <- if (basis == "c") chebMap(x, l, u) else x
  calcFunc <- switch(EXPR = basis, m = polyCalc, chebCalc)
  ret <- calcFunc(z, R$a)
  if ("b" %in% names(R)) {
    ret <- ret / calcFunc(z, R$b)
  }

  ret
}

# Function to calculate error between known and calculated values.
remErr <- function(x, R, fn, relErr, basis, l, u) {
  if (relErr) {
    y <- callFun(fn, x)
    (evalFunc(x, R, basis, l, u) - y) / y
  } else {
    evalFunc(x, R, basis, l, u) - callFun(fn, x)
  }
}

# Function to identify roots of the error equation for use as bounds in finding
# the maxima and minima.
findRoots <- function(x, R, fn, relErr, basis, l, u) {
  r <- double(length(x) - 1L)
  for (i in seq_along(r)) {
    intv <- c(x[i], x[i + 1L])
    root <- tryCatch(uniroot(remErr, interval = intv, extendInt = "no", R = R,
                             fn = fn, relErr = relErr, basis = basis,
                             l = l, u = u, tol = sqrt(.Machine$double.eps)),
                     error = function(cond) simpleError(trimws(cond$message)))

    # If there is no root in the interval, take the endpoint closest to zero.
    if (inherits(root, "simpleError")) {
      r[i] <- intv[which.min(abs(intv))]
    } else {
      r[i] <- root$root
    }
  }

  r
}

# F9 fix (hardening): compute the ZeroBasis perturbation for a candidate
# extremum x_i that landed exactly on a zero of fn (relative error is
# undefined there). The original fixed 1e-12 absolute nudge is a genuine
# no-op once |x_i| >~ 4.5e3 (a double's ulp there already exceeds 1e-12),
# silently leaving the zero-division problem in place -- at endpoints
# (x_i + 1e-12 == x_i) and in the interior alike (both candidates collapse
# to the same value). The fix escalates to a magnitude-scaled step ONLY
# when the plain absolute step would not move x_i; at ordinary |x_i| the
# step is byte-identical to the released 1e-12 absolute behavior (this
# matters: the ZeroBasis end-to-end cases run at the machine-precision
# floor, where changing the step's magnitude or direction perturbs
# convergence). x_i == 0 keeps the pure absolute step. Extracted from
# switchX so it is unit-testable independent of the optimizer.
#
# stepInto: additive step of magnitude >= peturb in the given direction
# (+1 = increase x_i, -1 = decrease), escalating to abs(x_i)*peturb only if
# the base peturb is absorbed. Used for the two endpoints (always stepping
# inward, into the interval).
stepInto <- function(x_i, dir, peturb = 1e-12) {
  s <- peturb
  if (x_i + dir * s == x_i) s <- abs(x_i) * peturb
  x_i + dir * s
}

zeroBasisPerturb <- function(x_i, l, u, fn, maximize, peturb = 1e-12) {
  if (x_i == l) {
    # Step inward from the lower endpoint (toward the interior).
    stepInto(x_i, 1, peturb)
  } else if (x_i == u) {
    # Step inward from the upper endpoint (mirror image of the l case).
    stepInto(x_i, -1, peturb)
  } else {
    # Interior: offer both directions and let the max/min pick. Escalate the
    # step magnitude only if the base absolute step is absorbed (large |x_i|).
    s <- peturb
    if (x_i - s == x_i || x_i + s == x_i) s <- abs(x_i) * peturb
    xreplace <- c(x_i - s, x_i + s)
    fnreplace <- callFun(fn, xreplace)
    if (maximize) {
      xreplace[which.max(fnreplace)]
    } else {
      xreplace[which.min(fnreplace)]
    }
  }
}

# Function to identify new x positions. This algorithm uses the multi-switch
# paradigm, not the single switch.
switchX <- function(r, l, u, R, fn, relErr, basis) {
  bottoms <- c(l, r)
  tops <- c(r, u)
  x <- double(length(bottoms))
  attr(x, "ZeroBasis") <- FALSE
  maximize <- sign(remErr(l, R, fn, relErr, basis, l, u)) == 1
  for (i in seq_along(x)) {
    intv <- c(bottoms[i], tops[i])
    # Tighter tolerances than the default lead to issues (AA: 2024-01-31).
    extrma <- tryCatch(optimize(remErr, interval = intv, R = R, fn = fn,
                                relErr = relErr, basis = basis, l = l, u = u,
                                maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    # If no extremum then the take endpoint with "better" value depending if we
    # are maximizing or minimizing.
    if (inherits(extrma, "simpleError")) {
      endPtErr <- remErr(intv, R, fn, relErr, basis, l, u)
      if (maximize) {
        x[i] <- intv[which.max(endPtErr)]
      } else {
        x[i] <- intv[which.min(endPtErr)]
      }
    } else {
      x[i] <- extrma[[1L]]
    }

    # Test endpoints for max/min even if an extremum was found.
    p <- c(bottoms[i], x[i], tops[i])
    E <- remErr(p, R, fn, relErr, basis, l, u)

    if (maximize) {
      x[i] <- p[which.max(E)]
    } else {
      x[i] <- p[which.min(E)]
    }

    # Test for 0 value at function if relative error
    if (relErr && callFun(fn, x[i]) == 0) {
      attr(x, "ZeroBasis") <- TRUE
      x[i] <- zeroBasisPerturb(x[i], l, u, fn, maximize)
    }

    # Flip maximize.
    maximize <- !maximize
  }

  x
}

# Check Remez iterations for convergence.
isConverged <- function(errs, expe, convrat, tol) {
  aerrs <- abs(errs)
  mxae <- max(aerrs)
  mnae <- min(aerrs)
  a_mxa_exp <- abs(mxae - expe)
  mx_mn <- mxae - mnae

  # Check observed errors are close enough to expected by ratio or tolerance.
  errDistance <- mxae / expe <= convrat ||
    (a_mxa_exp <= tol && a_mxa_exp > .Machine$double.eps)

  # Check observed errors are close enough to each other by ratio or tolerance.
  errMagnitude <- mxae / mnae <= convrat ||
    (mx_mn <= tol && mx_mn > .Machine$double.eps)

  # Converged if magnitude and distance are close and error oscillates in sign.
  isOscil(errs) && errDistance && errMagnitude
}

isUnchanging <- function(errs, errs_last, convrat, tol) {
  denomProblem <- which(errs_last == 0)
  # If any are actually 0, then perturb them by 1e-12. Ratio becomes 1 and
  # difference remains 0.
  if (length(denomProblem) > 0L) {
    errs[denomProblem] <- errs[denomProblem] + 1e-12
    errs_last[denomProblem] <- errs_last[denomProblem] + 1e-12
  }
  errsDiff <- abs(errs - errs_last)
  ratio <- abs(errs / errs_last)
  # F7 fix: the ratio test was one-sided (<= convrat), so errors shrinking
  # rapidly (e.g. 10x per iteration) satisfied it and were flagged as
  # "unchanging" -- premature stop on genuine improvement. The test must be
  # two-sided: only a ratio close to 1 (in EITHER direction) indicates
  # stagnation. The zero-denominator perturbation and absolute-difference
  # clause are unchanged.
  all(ratio <= convrat & ratio >= 1 / convrat) ||
    (all(errsDiff <= tol) && all(errsDiff > .Machine$double.eps))
}

# Check denominator polynomial for zero in the requested range.
checkDenom <- function(a, l, u, basis) {
  # M6 (F5 Option A): uniroot searches over raw x in [l, u]; wrapping the
  # evaluation function to map x -> z internally means uniroot's returned
  # root is ALREADY in raw x -- no back-conversion needed, and the error
  # message in remRat (which reports this root directly) stays correct
  # without any change there.
  calcFn <- if (basis == "m") {
    polyCalc
  } else {
    function(x, a) chebCalc(chebMap(x, l, u), a)
  }
  dngrRt <- tryCatch(uniroot(calcFn, c(l, u), extendInt = "no", a = a,
                             tol = .Machine$double.eps),
                     error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(dngrRt, "simpleError")) {
    return(NULL)
  } else {
    return(dngrRt$root)
  }
}

# Basis-aware per-coefficient scale bounding a_k's contribution to the
# approximation, for k = 0 .. n-1 (n = length(a)). Used by checkIrrelevant
# (F6) and tailContribution (F3).
# Monomial: unchanged xmax^k (xmax = max(|l|, |u|)).
# M6 (F5 Option A): Chebyshev basis now ALWAYS evaluates T_k on the mapped
# [-1, 1] domain, so |T_k(z)| <= 1 unconditionally -- the bound on a_k's
# contribution collapses to the trivial 1 for every k, exactly as flagged by
# the M3 comment this replaces. l/u are accepted but unused in this branch
# (kept for signature symmetry with the monomial branch and existing callers).
basisScale <- function(n, l, u, basis) {
  if (basis == "m") {
    xmax <- max(abs(l), abs(u))
    xmax ^ (seq_len(n) - 1L)
  } else {
    rep(1, n)
  }
}

# Check for coefficient irrelevancy.
checkIrrelevant <- function(a, l, u, zt, basis) {
  n <- length(a)
  if (!is.null(zt) && n > 0) {
    a <- ifelse(abs(a * basisScale(n, l, u, basis)) <= zt, 0, a)
  }

  a
}

# F3 fix: the n+1-restart "effectively zero" test in minimaxApprox() checked
# (a_n * xmax^(n-1)) > tailtol with no abs() on a_n, so ANY negative top
# coefficient -- of arbitrary magnitude -- passed and was silently dropped.
# Extracted as its own function (unit-testable) taking abs() of the
# coefficient and the basis-correct scale from basisScale (monomial:
# xmax^(n-1); Chebyshev: the F6 endpoint/unit-extremum bound of |T_{n-1}|).
tailContribution <- function(a_n, n, l, u, basis) {
  abs(a_n) * basisScale(n, l, u, basis)[n]
}
