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
  # M5 (barycentric): basis "b" evaluates the trial/fitted polynomial directly
  # from its stored barycentric representation R$bary = list(x, w, p) via the
  # second barycentric formula (baryEval). This is the accurate evaluation path
  # and the one findRoots/switchX/remErr/plot/minimaxErr all reach through this
  # single dispatch site, so those functions need no barycentric-specific code.
  # The field is named "bary" (not "b") so it does not collide with the
  # rational-denominator detection `"b" %in% names(R)` below.
  if (basis == "b") {
    return(baryEval(x, R$bary$x, R$bary$w, R$bary$p))
  }
  # M6 (F5 Option A): l/u are now REQUIRED (not defaulted) on purpose.
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

# CP-1 (E5 Phase 1): reference-local certificate threshold, shared by all four
# approximation paths. Calibrated on the rational-barycentric path (M5P2 F.2:
# healthy converged fits <= 1 + 1.2e-5, known sub-optimal fixed points
# >= 1.11) and re-verified on the polynomial paths (MP2 B1-1: healthy
# <= 1 + 3.4e-5, bad cases >= 1.77): 1.001 sits in a four-decade dead zone.
REFLOCALTOL <- 1.001

# CP-1: certificate check run at a CONVERGED exit only. The Remez loop's
# convergence test evaluates the trial at its own reference, where the trial
# is leveled by construction (or leveled up to solve residual on the classical
# paths), so it carries no off-reference information -- any exchange stall is
# certified as convergence (MP2 B1-1 mechanism D3). This check supplies the
# missing certificate: compare the returned approximation's dense-grid sup
# error against its leveled error; a ratio above REFLOCALTOL means the
# leveled error is only a lower bound (de la Vallee Poussin) on the true
# minimax error, not a certificate, and the caller raises the refLocal
# warning centrally.
#
# Skips (returns gridSup = NA, refLocal = FALSE, i.e. no certificate rather
# than a spurious one):
# * relErr with an exact zero of fn on the probe grid: the relative error
#   diverges as 1/fn around the zero, so the grid sup is dominated by the
#   blow-up and any threshold fires spuriously. Excluding only the exact-zero
#   points does not help -- neighboring points still carry the divergence.
#   (Same guard and rationale as interpRescue, M4.)
# * Floor gate: expe at or below ~100*eps*max(1, ||f||) (relErr: 100*eps).
#   Below that the "excess" of gridSup over expe is linear-solve residual
#   noise at magnitude ~O(10)*eps*||f||, not a basin miss (measured: the
#   exp[5, 6] degree-10 relaxed-convrat fit, expe 2.94e-12 vs gate 8.96e-12,
#   ratio 1.05 of pure noise). True sub-optimal fixed points in the floor
#   neighborhood sit well ABOVE the gate (sin degree 9, expe 3.2e-13 vs gate
#   2.2e-14, ratio 406) and still fire.
# * A non-finite probe evaluation (defensive; a converged fit has already
#   evaluated fn across the interval).
#
# The rational-barycentric path (remBaryRat) retains its own in-loop F.2
# check -- permanent there because non-normal rational problems have true
# alternants exceeding m + n + 2 points, which no exchange rule can hold --
# and shares only the REFLOCALTOL constant, not this gate.
refLocalCheck <- function(R, fn, relErr, basis, lower, upper, expe) {
  none <- list(gridSup = NA_real_, refLocal = FALSE)
  probe <- seq(lower, upper, length.out = 2001L)
  fv <- callFun(fn, probe)
  if (relErr && any(fv == 0)) return(none)
  normf <- max(abs(fv))
  if (!is.finite(normf)) return(none)                               # nocov
  gate <- 100 * .Machine$double.eps * if (relErr) 1 else max(1, normf)
  if (expe <= gate) return(none)
  e <- evalFunc(probe, R, basis, lower, upper) - fv
  if (relErr) e <- e / fv
  gridSup <- max(abs(e))
  if (!is.finite(gridSup)) return(none)                             # nocov
  # Absolute-excess condition (approved refinement to the ratio test): the
  # ratio detects the reference-local phenomenon, but when expe is only
  # modestly above the floor gate, a ratio a whisker over REFLOCALTOL can
  # correspond to an ABSOLUTE excess of a few eps*||f|| -- pure solve /
  # evaluation noise, platform-flappy by construction (measured: exp deg 11
  # Chebyshev, ratio 1.0011, excess 1.1e-15 ~ 4*eps*||f||). Require the
  # excess itself to be resolvable above floating-point noise. Every
  # ordinary-magnitude signal has excess many orders above this (B1-1
  # quartet: 1.3e-10 .. 3.1e-3), so no signal is lost.
  noiseFloor <- 10 * .Machine$double.eps * if (relErr) 1 else max(1, normf)
  list(gridSup = gridSup,
       refLocal = gridSup > REFLOCALTOL * expe &&
         gridSup - expe > noiseFloor)
}

# E5: identify ALL roots of the error equation over the whole interval
# [l, u]. The pre-E5 implementation searched only the length(x) - 1 intervals
# BETWEEN consecutive current reference points, so [l, x_1] and [x_last, u]
# were never searched (MP2 B1-1 mechanism D1): the root count was always
# exactly length(x) - 1, downstream switchX always returned exactly length(x)
# points, and a sign change between an endpoint and the outermost reference
# point was structurally invisible -- the root cause of the silent
# reference-local convergence family. Now the error is evaluated on an
# oversampled grid (10x the reference size, minimum 129 points -- a degree-n
# minimax error curve has at most ~n + 2 extrema, so 10x is generous), the
# current reference points are appended (their leveled +-h values anchor
# their sign runs at no extra cost), every strict sign change is refined by
# uniroot on its bracketing grid subinterval, and exact grid zeros are taken
# as roots directly. The root count is now variable -- that is the point;
# the old "no root -> endpoint closest to zero" substitution is gone (a
# degenerate curve simply yields fewer roots, handled by switchX's
# fallback). Signature and return type (sorted numeric vector) unchanged.
findRoots <- function(x, R, fn, relErr, basis, l, u) {
  nGrid <- max(129L, 10L * length(x))
  # Clamp the reference to [l, u] defensively (production references always
  # lie within; direct calls may not).
  g <- sort(unique(c(seq(l, u, length.out = nGrid), x[x >= l & x <= u])))
  e <- remErr(g, R, fn, relErr, basis, l, u)
  s <- sign(e)

  # Exact zeros on the grid are roots as-is (also covers tangential zeros,
  # which have no sign CHANGE to bracket).
  r <- g[s == 0]

  # Strict sign changes: refine each by uniroot on its grid bracket. The
  # bracket is verified to change sign, so uniroot cannot fail on its
  # f(lower)/f(upper) precondition; the tryCatch is defense in depth only.
  chg <- which(s[-length(s)] * s[-1L] < 0)
  for (i in chg) {
    intv <- c(g[i], g[i + 1L])
    root <- tryCatch(uniroot(remErr, interval = intv, extendInt = "no", R = R,
                             fn = fn, relErr = relErr, basis = basis,
                             l = l, u = u,
                             tol = sqrt(.Machine$double.eps))$root,
                     error = function(cond) {                       # nocov
                       intv[which.min(abs(e[c(i, i + 1L)]))]        # nocov
                     })                                             # nocov
    r <- c(r, root)
  }

  sort(unique(r))
}

# E5: shared window-selection rule for choosing N reference points from an
# ordered list of alternating extremum candidates. Among the size-N windows
# of consecutive candidates, choose the one maximizing the smallest |error|
# (FNT 2018 Step 3, "largest values satisfying the alternation"; ties break
# to the first window, preserving baryRatInitRef's pre-E5 which.max
# behavior bitwise). requireMax additionally restricts the windows to those
# CONTAINING the global |error| argmax -- the classical exchange invariant
# (the reference must hold the global extremum for h to increase strictly
# toward E_true). baryRatInitRef calls with requireMax = FALSE (its pre-E5
# semantics, an initialization heuristic where the invariant is not needed);
# switchX calls with requireMax = TRUE. Note the two differ materially:
# max-min windowing alone can DROP the global max (e.g. |e| =
# {10, 1e-6, 5, 5, 5, 5}, N = 4 selects the last four).
selectAlternantWindow <- function(vals, N, requireMax) {
  nw <- length(vals) - N + 1L
  wmin <- vapply(seq_len(nw), function(i) min(vals[i:(i + N - 1L)]), double(1L))
  if (requireMax) {
    imax <- which.max(vals)
    ok <- seq_len(nw) <= imax & seq_len(nw) + N - 1L >= imax
    wmin[!ok] <- -Inf
  }
  which.max(wmin)
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

# E5: identify new reference positions from the error curve. Multi-switch
# exchange, redesigned error-curve-first (PT09/Chebfun shape). The pre-E5
# implementation seeded a maximize/minimize schedule from sign(remErr(l)) and
# flipped it mechanically per bracket; when the error curve had one more
# oscillation than the bracket structure assumed, the schedule desynced and
# the exchange returned roots instead of extrema, destroying the reference
# around a correctly captured global max (MP2 B1-1 mechanism D2, traced).
# Now each bracket's direction comes from the ERROR'S OWN SIGN at the
# bracket midpoint, so it cannot desync, and the candidate count is variable
# (one per sign-region of the curve, endpoints included).
#
# New argument xk = the CURRENT reference (no default: any missed caller
# fails loudly). It supplies the target size N = length(xk) and the
# degenerate fallback: with fewer than N alternating candidates (floor
# regime, flat error curve) the previous reference is returned UNCHANGED, so
# the next trial reproduces itself and the loop exits through the existing
# isUnchanging stall machinery -- deliberately NO new degenerate handling
# (M3 lesson: never intercept or condition the singular/restart flow).
#
# Selection (the organ the pre-E5 pipeline lacked, since its candidate count
# was structurally fixed): consecutive same-sign candidates are merged
# keeping the larger |error| (tangential zeros can produce same-sign
# neighbors); a surplus is resolved by selectAlternantWindow with
# requireMax = TRUE, so the returned reference always alternates, always
# contains the global |error| extremum, and maximizes its smallest |error|
# -- the classical conditions under which the leveled error increases
# strictly toward the true minimax.
switchX <- function(r, l, u, R, fn, relErr, basis, xk) {
  N <- length(xk)
  brk <- sort(unique(c(l, r, u)))
  nb <- length(brk) - 1L
  cx <- ce <- double(nb)

  for (i in seq_len(nb)) {
    intv <- c(brk[i], brk[i + 1L])
    # Direction from the error's own sign in this bracket (midpoint sample;
    # the sign is constant between consecutive roots).
    maximize <- remErr((intv[1L] + intv[2L]) / 2, R, fn, relErr,
                       basis, l, u) > 0
    # Tighter tolerances than the default lead to issues (AA: 2024-01-31).
    extrma <- tryCatch(optimize(remErr, interval = intv, R = R, fn = fn,
                                relErr = relErr, basis = basis, l = l, u = u,
                                maximum = maximize),
                       error = function(cond) simpleError(trimws(cond$message)))

    # If no extremum then take the endpoint with "better" value depending if
    # we are maximizing or minimizing.
    if (inherits(extrma, "simpleError")) {
      endPtErr <- remErr(intv, R, fn, relErr, basis, l, u)            # nocov
      xi <- intv[if (maximize) which.max(endPtErr) else                # nocov
        which.min(endPtErr)]                                # nocov
    } else {
      xi <- extrma[[1L]]
    }

    # Test endpoints for max/min even if an extremum was found.
    p <- c(intv[1L], xi, intv[2L])
    E <- remErr(p, R, fn, relErr, basis, l, u)
    j <- if (maximize) which.max(E) else which.min(E)
    cx[i] <- p[j]
    ce[i] <- E[j]
  }

  # Endpoints are ALWAYS standalone candidates (not merely contestants in
  # their brackets' 3-point check). At a parity-degenerate trial (even/odd
  # fn on a symmetric reference, where the leveled h is legitimately ~0 and
  # the trial is essentially the interpolant), the error curve's interior
  # lobes offer only N - 1 alternating extrema -- the endpoint lobes are
  # pinned to zero by the endpoint nodes -- and the exchange would stall on
  # the previous reference forever (measured: cos deg 4 barycentric). The
  # endpoints' leveled values -sigma_i*h, though eps-scale, carry the
  # alternation signs the lobe list lacks: after the same-sign merge below,
  # exactly the endpoint whose sign OPPOSES its neighboring lobe survives,
  # restoring an N-point alternating candidate set (min |e| = h, global max
  # retained -- a legal exchange) whose asymmetry breaks the parity
  # degeneracy on the next trial.
  eEnd <- remErr(c(l, u), R, fn, relErr, basis, l, u)
  if (cx[1L] != l) {
    cx <- c(l, cx)
    ce <- c(eEnd[1L], ce)
    nb <- nb + 1L
  }
  if (cx[length(cx)] != u) {
    cx <- c(cx, u)
    ce <- c(ce, eEnd[2L])
    nb <- nb + 1L
  }

  # Merge consecutive same-sign candidates, keeping the larger |error|.
  keep <- logical(nb)
  keep[1L] <- TRUE
  last <- 1L
  for (i in seq_len(nb)[-1L]) {
    if (sign(ce[i]) != sign(ce[last])) {
      keep[i] <- TRUE
      last <- i
    } else if (abs(ce[i]) > abs(ce[last])) {
      keep[last] <- FALSE
      keep[i] <- TRUE
      last <- i
    }
  }
  cx <- cx[keep]
  ce <- ce[keep]

  # Degenerate fallback: too few alternating candidates to fill the
  # reference. Return the previous reference unchanged and let the existing
  # stall machinery exit the loop.
  if (length(cx) < N) {
    x <- xk
    attr(x, "ZeroBasis") <- FALSE
    return(x)
  }

  if (length(cx) > N) {
    i0 <- selectAlternantWindow(abs(ce), N, requireMax = TRUE)
    sel <- i0:(i0 + N - 1L)
    cx <- cx[sel]
    ce <- ce[sel]
  }

  x <- cx
  attr(x, "ZeroBasis") <- FALSE

  # Test for 0 value at function if relative error (unchanged contract).
  if (relErr) {
    for (i in seq_len(N)) {
      if (callFun(fn, x[i]) == 0) {
        attr(x, "ZeroBasis") <- TRUE
        x[i] <- zeroBasisPerturb(x[i], l, u, fn, ce[i] > 0)
      }
    }
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
