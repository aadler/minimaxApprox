# Copyright Avraham Adler (c) 2026
# SPDX-License-Identifier: MPL-2.0+

# Module M5, Phase 1: barycentric-Remez polynomial approximation (Pachon &
# Trefethen 2009, BIT 49:721-741), found at
# https://www.chebfun.org/publications/remez.pdf. This path solves NO linear
# system: the trial polynomial is represented by its values at the reference
# (barycentric Lagrange form) and the leveled error has a closed form, so the
# conditioning failures of the classical bases (F4 degree ceiling, F5 off-range
# accuracy) do not arise here. It therefore does NOT route through the classical
# singular/restart machinery or interpRescue (M4) -- see minimaxApprox() where
# basis == "b" is dispatched and the two rescue blocks are guarded off.

# Second-kind Chebyshev points (extrema; endpoints INCLUDED), sorted ascending.
# The package's chebNodes() are first-kind (roots; endpoints excluded); PT09
# initializes the reference at these N = degree + 2 extrema. Kept separate so
# the classical paths' release-baselined initial reference is untouched.
chebNodes2 <- function(n, a, b) {
  n <- as.integer(n)
  sort(0.5 * (a + b) + 0.5 * (b - a) * cos((seq_len(n) - 1L) * pi / (n - 1L)))
}

# Separate any coincident (or sub-tolerance) reference nodes. switchX can return
# two reference points at the same x. For an even/odd function, the symmetric
# initial reference makes the closed-form leveled error cancel to machine noise
# (mathematically h = 0 there, though fn is NOT in P_n), so the trial polynomial
# nearly interpolates fn at the nodes. Therefore, the error curve then grazes
# zero AT a node and findRoots brackets a spurious root there, collapsing two
# reference points. The classical QR solve tolerates a near-duplicate row;
# baryWeights cannot, since a zero gap gives log(0) = -Inf -> Inf weight -> NaN
# leveled error. Nudging the later of each coincident pair apart by a tiny
# fraction of the interval breaks the symmetry: the next iteration's leveled
# error becomes meaningful and normal convergence resumes (verified on
# tanh(x+.5)-tanh(x-.5) deg 10, which reaches the PT09 Table-1 oracle
# 3.0009195e-7). x is sorted ascending (switchX order). Genuine Remez extrema
# are never this close (minimum spacing ~ (u-l)/n^2), so the guard only ever
# fires on the degenerate collapse, never on real points.
separateNodes <- function(x, l, u) {
  sep <- 1e-10 * (u - l)
  for (k in seq_len(length(x) - 1L)) {
    if (x[k + 1L] - x[k] < sep) {
      x[k + 1L] <- min(x[k] + sep, u)
    }
  }
  # A collision at (or pushed to) the UPPER endpoint cannot be separated by
  # the upward pass: the min(, u) clamp pins the later point at u and leaves
  # the pair coincident. Sweep downward from the top for any remaining
  # collision, pushing the EARLIER point down instead. For every input the
  # forward pass already separates (all healthy references and all interior
  # collapses), this loop compares and moves nothing, so previously-working
  # behavior is bitwise unchanged; it acts only where the old code returned
  # coincident nodes, which downstream turned into Inf barycentric weights
  # (polynomial path) or a non-finite basis matrix (rational path).
  for (k in rev(seq_len(length(x) - 1L))) {
    if (x[k + 1L] - x[k] < sep) {
      x[k] <- max(x[k + 1L] - sep, l)
    }
  }
  x
}

# Capacity-scaled barycentric weights over ALL points of x (paper 3.4):
# w_j = prod_{v!=j} sign(x_j - x_v) / exp(sum_{v!=j} log|x_j - x_v| / C) with
# capacity C = (u - l) / 4 dividing every difference. Assembled as a sign
# product times exp(-sum(log|.|)) so wide intervals and high degree do not
# over/underflow. Because the leveled error (levelError) and the barycentric
# evaluation (baryEval) are each RATIOS of the weights, any common scale factor,
# including the capacity scaling itself, cancels exactly. The scaling is purely
# for floating-point range, never for the mathematical result.
baryWeights <- function(x, l, u) {
  cap <- (u - l) / 4
  n <- length(x)
  w <- double(n)
  for (j in seq_len(n)) {
    d <- (x[j] - x[-j]) / cap
    w[j] <- prod(sign(d)) * exp(-sum(log(abs(d))))
  }
  w
}

# Closed-form leveled error (paper eq. 3.9 for absolute error):
#   h = sum(w f) / sum(sigma w).
# Relative-error form, derived from the relative equioscillation condition
# (f_i - p_i) / f_i = sigma_i h  (=> p_i = f_i (1 - sigma_i h)) together with
# the degree-n condition sum_j w_j p(x_j) = 0:
#   sum_j w_j f_j (1 - sigma_j h) = 0  =>  h = sum(w f) / sum(sigma w f).
# (Verified numerically by equioscillation of the RELATIVE error itself, not by
# comparison to the absolute optimum -- the two optima differ.)
levelError <- function(w, sigmaB, f, relErr) {
  if (relErr) {
    sum(w * f) / sum(sigmaB * w * f)
  } else {
    sum(w * f) / sum(sigmaB * w)
  }
}

# Second barycentric formula (paper 3.4) evaluating the trial polynomial p
# (represented by node values pval at nodes xref with weights w) at arbitrary
# points x. The formula is 0/0 exactly at a node, so each requested point that
# coincides with a reference node short-circuits to that node's value.
baryEval <- function(x, xref, w, pval) {
  vapply(x, function(xk) {
    d <- xk - xref
    hit <- which(d == 0)
    if (length(hit) > 0L) {
      pval[hit[1L]]
    } else {
      wd <- w / d
      sum(wd * pval) / sum(wd)
    }
  }, double(1L))
}

# Build the trial polynomial (in barycentric form) and leveled error for a
# reference x. Returns the R-object shape consumed by evalFunc's basis == "b"
# branch: R = list(bary = list(x, w, p)). The `bary` name (not `b`) is
# deliberate -- evalFunc detects the rational denominator via `"b" %in%
# names(R)`, so a field literally named `b` would collide.
baryTrial <- function(x, fn, relErr, l, u, sigmaB) {
  x <- separateNodes(x, l, u)
  w <- baryWeights(x, l, u)
  f <- callFun(fn, x)
  zeroBasis <- FALSE

  # Leveled error via the shared closed form (levelError is the single source
  # of truth for both the absolute and relative branches; baryTrial no longer
  # duplicates the formula).
  h <- levelError(w, sigmaB, f, relErr)

  # relErr guard: a near-zero denominator (sum(sigma * w * f), the quantity
  # levelError divides by in relErr mode) means fn has (near) a zero close to
  # the reference, where relative error is ill-defined -- the same pathology
  # the classical ZeroBasis machinery (switchX/zeroBasisPerturb) handles. Flag
  # it (surfaces as the existing ZeroBasis warning in minimaxApprox()) and keep
  # h finite so the iteration can still report a result rather than crashing.
  if (relErr) {
    d <- sigmaB * w * f
    if (!is.finite(h) || abs(sum(d)) <= .Machine$double.eps * max(abs(d))) {
      zeroBasis <- TRUE
      # NaN reach the `h == 0` test below and abort the whole approximation. In
      # relErr mode the guard above already zeroes a non-finite h; in absolute
      # mode sum(sigma * w) does not vanish for a distinct-node reference
      # (separateNodes guarantees distinctness), so this is belt-and-suspenders
      # with no constructible trigger -- excluded from coverage.
      if (!is.finite(h)) h <- 0                                  # nocov

    }
  }

  # Absolute-path defense in depth: separateNodes removes the coincident-node
  # cause of a non-finite h, but should the leveled error still come back
  # non-finite for any reason, fall back to 0 (interpolant) rather than let a
  # NaN reach the `h == 0` test below and abort the whole approximation.
  if (!is.finite(h)) h <- 0                                     # nocov

  # paper 3.6: if h is exactly 0, perturb it so the trial polynomial is not
  # identically f at the reference (which would make the next exchange stall).
  # For a function already resolved to precision, remBary short-circuits on the
  # h-floor BEFORE reaching a loop that needs this, so this only matters mid-
  # iteration.
  if (h == 0) h <- 1e-19

  p <- if (relErr) f * (1 - sigmaB * h) else f - sigmaB * h
  list(R = list(bary = list(x = x, w = w, p = p)), h = h, x = x,
       zeroBasis = zeroBasis)
}

# One-point (first Remez algorithm) exchange, paper 2.2, used only as the
# overshoot safeguard: swap a single reference point for the global extremum of
# the current error curve, preserving sign alternation. Reuses findRoots to
# bracket the error's critical points, then locates the extremum in each
# bracket. Rare path (triggered when the full multi-switch exchange picks a
# large-norm trial polynomial), so clarity is favored over micro-optimization.
onePointExchange <- function(xk, R, fn, relErr, l, u) {
  err <- function(z) remErr(z, R, fn, relErr, "b", l, u)
  r <- findRoots(xk, R, fn, relErr, "b", l, u)
  brk <- sort(unique(c(l, r, u)))
  ext <- double(0L)
  for (i in seq_len(length(brk) - 1L)) {
    br <- c(brk[i], brk[i + 1L])
    if (br[2L] <= br[1L]) next             # nocov Cannot occur.
    m <- tryCatch(optimize(function(z) abs(err(z)), br, maximum = TRUE)$maximum,
                  error = function(e) br[which.max(abs(err(br)))])
    cand <- c(br, m)
    ext <- c(ext, cand[which.max(abs(err(cand)))])
  }
  ext <- sort(unique(ext))
  eext <- err(ext)
  gi <- which.max(abs(eext))
  xnew <- ext[gi]
  snew <- sign(eext[gi])
  exk <- err(xk)
  sxk <- sign(exk)
  nk <- length(xk)

  # Choose the point to drop so the alternation is maintained (paper 2.2).
  if (xnew > xk[1L] && xnew < xk[nk]) {
    same <- which(sxk == snew)
    xold <- xk[same[which.min(abs(xk[same] - xnew))]]
  } else if (xnew <= xk[1L]) {
    xold <- if (sxk[1L] == snew) xk[1L] else xk[nk]
  } else {
    xold <- if (sxk[nk] == snew) xk[nk] else xk[1L]
  }

  sort(c(xnew, xk[xk != xold]))
}

# Assemble the completed barycentric result. Beyond the classical fields (a,
# expe, mxae, i, x, converged, ...), it recovers the familiar coefficient
# vectors for continuity with the other bases and reports the conversion
# residual:
#  - a: Chebyshev coefficients in the MAPPED variable (M6 semantics), obtained
#    by evaluating the (stable) barycentric polynomial at n + 1 first-kind
#    Chebyshev nodes and solving the well-conditioned mapped Chebyshev-
#    Vandermonde system. aMono is then produced downstream in minimaxApprox()
#    by the existing cheb2mon + composeAffine pipeline, exactly as for basis
#    "c".
#  - convResid: max |p_bary - p_viaCoeffs| on a 1001-point probe grid. The
#    barycentric evaluation is the accurate representation; the coefficient
#    conversion is the one step that can lose accuracy at high degree, so its
#    gap is measured and stored (documented caveat).
finishBary <- function(trial, x, expe, mxae, i, converged, unchanged,
                       unchanging_i, zeroBasisError, l, u, n) {
  bx <- trial$R$bary$x
  bw <- trial$R$bary$w
  bp <- trial$R$bary$p

  xa <- chebNodes(n + 1L, l, u)
  va <- baryEval(xa, bx, bw, bp)
  a <- solve(chebMat(chebMap(xa, l, u), n), va)

  g <- seq(l, u, length.out = 1001L)
  convResid <- max(abs(baryEval(g, bx, bw, bp) - chebCalc(chebMap(g, l, u), a)))

  list(a = a, bary = list(x = bx, w = bw, p = bp), convResid = convResid,
       expe = expe, mxae = mxae, i = i, x = x, converged = converged,
       unchanged = unchanged, unchanging_i = unchanging_i,
       zeroBasisError = zeroBasisError)
}

# Main barycentric-Remez driver (polynomial). Mirrors remPoly's loop structure
# so the opts contract (maxiter, miniter, conviter, convrat, tol, showProgress)
# and the shared isConverged/isUnchanging behave identically across bases.
# tailtol/ztol are classical-path concepts (no coefficient vector to zero mid-
# iteration, no singular restart) and are simply not consulted here.
remBary <- function(fn, lower, upper, degree, relErr, opts) {
  n <- as.integer(degree)
  N <- n + 2L
  sigmaB <- (-1) ^ (seq_len(N) - 1L)
  relErrZeroBasis <- FALSE

  # ||f|| estimate (for the overshoot safeguard) and the machine-precision
  # floor used to detect the F4 family (function resolved to precision at the
  # requested degree).
  probe <- seq(lower, upper, length.out = 2001L)
  normf <- max(abs(callFun(fn, probe)))
  if (!is.finite(normf) || normf == 0) normf <- 1
  hFloor <- 10 * .Machine$double.eps * if (relErr) 1 else max(1, normf)

  # Initial reference: second-kind Chebyshev points.
  x <- chebNodes2(N, lower, upper)
  trial <- baryTrial(x, fn, relErr, lower, upper, sigmaB)
  x <- trial$x
  relErrZeroBasis <- relErrZeroBasis || trial$zeroBasis

  errs_last <- remErr(x, trial$R, fn, relErr, "b", lower, upper)
  mxae <- max(abs(errs_last))

  # relErr node-on-exact-zero guard (M5). When relErr is requested and fn is
  # EXACTLY 0 (to the bit) at a reference node, the trial value there is
  # p = f * (1 - sigma h) = 0 as well, so the relative error is (0 - 0)/0 = NaN.
  # The forced second-kind endpoint nodes make this the common case whenever an
  # endpoint is an exact zero of fn (e.g. sin(0) on [0, b]). A NaN error curve
  # cannot be exchanged: switchX errors with "missing value where TRUE/FALSE
  # needed". This is distinct from the WELL-POSED relErr case where fn merely
  # has a zero NEAR (or at an interior node that is only floating-point-near) a
  # reference -- there f and p are tiny-but-nonzero, the ratio is finite or Inf,
  # and the iteration tolerates it and converges (verified: cos[0,pi] with a
  # node at pi/2 ~ 6.12e-17, sin[-1,1] with a node at 0 ~ 6.12e-17). Detect the
  # true 0/0 case SPECIFICALLY via is.nan (NOT is.finite, which would also fire
  # on the tolerable Inf) and stop with a clear message: the relative-error
  # minimax genuinely does not exist at a node-on-exact-zero.
  if (relErr && anyNA(errs_last) && any(is.nan(errs_last))) {
    stop("Relative error is undefined because 'fn' is exactly 0 at a ",
         "reference  point (typically an endpoint of the [lower, upper] ",
         "interval). Use  absolute error (relErr = FALSE), or choose an ",
         "interval whose endpoints  are not zeros of 'fn'.", call. = FALSE)
  }

  expe <- abs(trial$h)
  converged <- unchanged <- FALSE
  unchanging_i <- i <- 0L

  # F4 family, natural in-basis fix (paper 3.6): if fn is (numerically) a
  # polynomial of degree <= n, the trial polynomial already interpolates it
  # and there is no work to do -- return E = 0 immediately, no iteration, no
  # restart, no interpRescue. The near-eps warning is raised centrally in
  # minimaxApprox().
  #
  # The detection MUST use the error on the fine probe grid, NOT the error at
  # the reference points. remErr at the reference equals -sigma_i * h by
  # construction, so max|err| there is just |h| = expe and carries no
  # independent information. For an even/odd fn on the symmetric initial
  # reference, sum(w*f) cancels and h is spuriously ~0 (measured 2.8e-17 for
  # tanh(x+.5)-tanh(x-.5), deg 10) even though the true deg-10 error is 7e-7.
  # The off-reference grid error is the only signal that separates a genuine
  # polynomial from an accidental symmetry cancellation.
  mxaeGrid <- max(abs(remErr(probe, trial$R, fn, relErr, "b", lower, upper)))
  if (expe <= hFloor && mxaeGrid <= hFloor) {
    return(finishBary(trial, x, expe, mxae, 0L, TRUE, FALSE, 0L,
                      relErrZeroBasis, lower, upper, n))
  }

  repeat {
    if (i >= opts$maxiter) break
    i <- i + 1L

    r <- findRoots(x, trial$R, fn, relErr, "b", lower, upper)
    xFull <- switchX(r, lower, upper, trial$R, fn, relErr, "b")
    zb <- attr(xFull, "ZeroBasis")

    # PT09 overshoot safeguard (paper 3.5): if the current trial polynomial's
    # error at the full-exchange reference blows up relative to ||f|| (a large-
    # norm trial polynomial, usually near an endpoint), discard the full
    # exchange and redo it from the PREVIOUS reference with a one-point swap.
    # Measured note: with this package's mature findRoots/switchX exchange, the
    # full-exchange step does not overshoot on standard inputs (max observed
    # ratio ~0.36 on Runge deg 10/20, |x| deg 11, tanh(20x) deg 16, etc. -- all
    # far below the 1e5 trigger). The master-pass prototype hit overshoot only
    # because it used a crude exchange. The branch is therefore defensive and
    # not reachable through the public API here; onePointExchange itself is
    # exercised by a direct unit test. Kept per PT09 for robustness.
    errFull <- remErr(xFull, trial$R, fn, relErr, "b", lower, upper)
    scale_b <- if (relErr) 1 else normf
    if (max(abs(errFull)) / scale_b > 1e5) {
      x <- onePointExchange(x, trial$R, fn, relErr, lower, upper) # nocov
    } else {
      x <- xFull
    }
    relErrZeroBasis <- relErrZeroBasis || isTRUE(zb)

    trial <- baryTrial(x, fn, relErr, lower, upper, sigmaB)
    x <- trial$x
    relErrZeroBasis <- relErrZeroBasis || trial$zeroBasis
    errs <- remErr(x, trial$R, fn, relErr, "b", lower, upper)
    mxae <- max(abs(errs))
    expe <- abs(trial$h)

    if (opts$showProgress) {
      message("i: ", i, " E: ", fC(expe), " maxErr: ", fC(mxae))
    }

    if (isConverged(errs, expe, opts$convrat, opts$tol) && i >= opts$miniter) {
      converged <- TRUE
      break
    }

    if (isUnchanging(errs, errs_last, opts$convrat, opts$tol)) {
      unchanging_i <- unchanging_i + 1L
      if (unchanging_i >= opts$conviter) {
        unchanged <- TRUE
        break
      }
    }

    errs_last <- errs
  }

  finishBary(trial, x, expe, mxae, i, converged, unchanged, unchanging_i,
             relErrZeroBasis, lower, upper, n)
}
