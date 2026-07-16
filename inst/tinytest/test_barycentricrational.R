# Copyright Avraham Adler (c) 2026
# SPDX-License-Identifier: MPL-2.0+

# M5 Phase 2: rational barycentric-Remez (FNT 2018 subset).

tol <- sqrt(.Machine$double.eps)
sW <- function(x) suppressWarnings(x)
sM <- function(x) suppressMessages(x)

nS <- getNamespace("minimaxApprox")
baryRatCtilde <- get("baryRatCtilde", nS)
baryRatSolve <- get("baryRatSolve", nS)
baryRatPQ <- get("baryRatPQ", nS)
baryRatPoleCheck <- get("baryRatPoleCheck", nS)
chebNodes2 <- get("chebNodes2", nS)
callFun <- get("callFun", nS)
baryEval <- get("baryEval", nS)
chebCalc <- get("chebCalc", nS)
chebMap <- get("chebMap", nS)

# ---- Structural identities (FNT Lemma 4 / Corollary 6) ---------------------
# The QR factor of the C-tilde matrix must be orthonormal AND S-orthogonal
# (Q1' S Q1 = 0, S = diag(+1, -1, ...)); the latter is the load-bearing
# identity that reduces the leveled-error computation to a SYMMETRIC
# eigenproblem. Both must hold at machine precision, and they pin the
# support-row sign convention in baryRatCtilde (a wrong sign leaves Q1
# orthonormal but visibly breaks nothing except downstream selection, so the
# identity is tested directly here).
x6 <- chebNodes2(6L, -1, 1)
t3 <- x6[c(2L, 4L, 6L)]
Q1 <- qr.Q(qr(baryRatCtilde(x6, t3, 0.5)))
sig <- (-1) ^ (0:5)
expect_true(max(abs(crossprod(Q1) - diag(3L))) < 1e-13)
expect_true(max(abs(crossprod(Q1, sig * Q1))) < 1e-13)

# The solve on a Chebyshev reference must satisfy the leveled-error
# equioscillation conditions AT the reference to machine precision:
# |r(x_i) - f(x_i)| = |h| for every reference point, with alternating signs.
sol <- baryRatSolve(x6, callFun(exp, x6), 2L, 2L, -1, 1)
expect_null(sol$fail)
r6 <- baryEval(x6, sol$t, sol$beta, sol$p)
expect_true(max(abs(abs(r6 - exp(x6)) - abs(sol$h))) < 1e-13)
expect_true(all(diff(sign(r6 - exp(x6))) != 0))

# ---- Oracle parity with the classical Cody path (healthy smooth cases) ----
# Cross-path tolerances are kept loose enough (1e-6 relative) to be robust
# to BLAS/LAPACK differences; both paths converge much closer than this on
# the reviewer platform (measured 1.4e-11 for (2, 2)).
b22 <- minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b")
c22 <- minimaxApprox(exp, -1, 1, c(2L, 2L))
expect_true(abs(b22$ExpErr - c22$ExpErr) / c22$ExpErr < 1e-6)
expect_false(b22$Warning)

b44 <- minimaxApprox(exp, -1, 1, c(4L, 4L), basis = "b")
c44 <- minimaxApprox(exp, -1, 1, c(4L, 4L))
expect_true(abs(b44$ExpErr - c44$ExpErr) / c44$ExpErr < 1e-4)

# Nondiagonal both ways. These also pin the Krylov ones-seed adjudication:
# with the paper's literal f-values seed the degree constraint silently
# fails and these parities are violated by orders of magnitude.
b32 <- minimaxApprox(exp, -1, 1, c(3L, 2L), basis = "b")
c32 <- minimaxApprox(exp, -1, 1, c(3L, 2L))
expect_true(abs(b32$ExpErr - c32$ExpErr) / c32$ExpErr < 1e-4)

b23 <- minimaxApprox(exp, -1, 1, c(2L, 3L), basis = "b")
c23 <- sW(minimaxApprox(exp, -1, 1, c(2L, 3L)))
expect_true(abs(b23$ExpErr - c23$ExpErr) / c23$ExpErr < 1e-4)

# Requested degrees are honored exactly in the recovered coefficients.
expect_length(b32$a, 4L)
expect_length(b32$b, 3L)
expect_length(b23$a, 3L)
expect_length(b23$b, 4L)
expect_equal(b23$b[1L], 1, tolerance = tol)

# ---- |x|: the case class the barycentric basis exists for ------------------
# Two independent algorithms (classical Cody and barycentric FNT) agree on
# E_[2,2](|x|), and the barycentric result equioscillates globally (checked
# on a dense grid, not just at the reference).
bax <- minimaxApprox(function(x) abs(x), -1, 1, c(2L, 2L), basis = "b")
expect_equal(bax$ExpErr, 4.36890126e-2, tolerance = 1e-7)
expect_true(bax$ObsErr / bax$ExpErr - 1 < 1e-10)
g <- seq(-1, 1, length.out = 20001L)
expect_true(max(abs(minimaxEval(g, bax) - abs(g))) / bax$ExpErr < 1 + 1e-6)

# ---- Exact representability (F4 family) ------------------------------------
# f exactly of the requested type: the leveled error sits at the machine-
# precision floor, the short-circuit returns the interpolant with zero
# iterations, and the standard near-eps warning fires.
bex <- sW(minimaxApprox(function(x) 1 / (1 + x ^ 2), -1, 1, c(0L, 2L),
                        basis = "b"))
expect_true(bex$ExpErr < 1e-13)
expect_identical(bex$iterations, 0L)
expect_true(bex$Warning)

# ---- Evaluation seams -------------------------------------------------------
# The stored bary element evaluates the SAME rational as the recovered
# mapped-Chebyshev a/b quotient and the monomial aMono/bMono quotient, and
# minimaxErr is bounded by the certified error on interior points.
z <- seq(-1, 1, length.out = 41L)
vb <- minimaxEval(z, b22)
expect_equal(vb, sM(minimaxEval(z, b22, basis = "c")), tolerance = 1e-12)
expect_equal(vb, sM(minimaxEval(z, b22, basis = "m")), tolerance = 1e-12)
expect_true(max(abs(minimaxErr(z, b22))) <= b22$ExpErr * (1 + 1e-9))
expect_true(all(names(b22$convResid) == c("a", "b")))
expect_true(all(b22$convResid < 1e-12))
expect_true(all(c("x", "w", "p") %in% names(b22$bary)))

# baryRatPQ residue formula: an exact support-point hit must agree with the
# limit from nearby points.
pq0 <- baryRatPQ(sol$t[2L], sol$t, sol$alpha, sol$beta, 0.5)
pq1 <- baryRatPQ(sol$t[2L] + 1e-9, sol$t, sol$alpha, sol$beta, 0.5)
expect_equal(pq0$q, pq1$q, tolerance = 1e-6)

# A one-signed denominator passes the pole check; an alternating-weight sign
# flip forced by hand is caught.
expect_null(baryRatPoleCheck(sol$t, sol$alpha, sol$beta, -1, 1))
expect_true(is.numeric(baryRatPoleCheck(sol$t, sol$alpha, abs(sol$beta),
                                        -1, 1)))

# ---- Input contracts --------------------------------------------------------
expect_error(minimaxApprox(exp, -1, 1, c(2L, 2L), relErr = TRUE,
                           basis = "b"),
             "Relative error is not yet supported for rational")
expect_error(minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b",
                           xi = c(-1, 0, 1)),
             "needs to have 6 elements")
bxi <- minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b",
                     xi = chebNodes2(6L, -1, 1))
expect_true(abs(bxi$ExpErr - c22$ExpErr) / c22$ExpErr < 1e-6)

# ---- Wide interval: capacity scaling ---------------------------------------
bw <- minimaxApprox(exp, 0, 10, c(3L, 2L), basis = "b")
expect_true(is.finite(bw$ExpErr))
expect_true(bw$ExpErr > 0)
gw <- seq(0, 10, length.out = 5001L)
expect_true(max(abs(minimaxEval(gw, bw) - exp(gw))) / bw$ExpErr < 1 + 1e-6)

# ---- Coverage-completion tests (post-covr review) ---------------------------
baryRatFail <- get("baryRatFail", nS)
baryRatInitRef <- get("baryRatInitRef", nS)
separateNodes <- get("separateNodes", nS)
fEC <- function(x) exp(cos(x))

# baryRatPoleCheck: a denominator weight of exactly 0 at a support point that
# is also a probe-grid point makes q hit 0.0 exactly (the s == 0 branch).
expect_equal(baryRatPoleCheck(c(-0.5, 0, 0.5), c(1, 1, 1), c(1, 0, 1),
                                  -1, 1), 0L, tolerance = tol)

# Coincident reference points must return the clean "collapse" failure from
# baryRatSolve (defense in depth), never a raw qr()/low-level error.
xdup <- c(-1, -0.5, -0.5, 0, 0.5, 1)
sdup <- baryRatSolve(xdup, callFun(exp, xdup), 2L, 2L, -1, 1)
expect_identical(sdup$fail, "collapse")

# separateNodes: a collision at the upper endpoint cannot be pushed upward
# (the clamp pins it at u); the backward pass must separate it downward.
xsep <- separateNodes(c(-1, 0, 1, 1), -1, 1)
expect_true(min(diff(xsep)) > 0)
expect_true(all(xsep >= -1 & xsep <= 1))

# All four detect-and-stop message arms, directly.
expect_error(baryRatFail("rank", 3L, 2L), "rank-deficient basis matrix")
expect_error(baryRatFail("nopolefree", 3L, 2L),
             "pole inside the approximation interval")
expect_error(baryRatFail("zeroweight", 3L, 2L),
             "denominator weight collapsed to zero")
expect_error(baryRatFail("collapse", 3L, 2L),
             "reference points collapsed onto each other")
expect_error(baryRatFail("rank", 3L, 2L), "c(2, 1)", fixed = TRUE)  # sugg degs

# Defect family (even function): detect-and-stop end-to-end. exp(cos(x)) is
# even, so its best rational approximations are non-normal at many requested
# types; the eigenpair selection correctly finds no pole-free solution. The
# (3, 3) case also exercises the AAA error-extrema window initialization
# (verified below to differ from the Chebyshev fallback); the (8, 8) case
# fails at a MID-ITERATION reference, exercising the in-loop dispatch.
expect_false(isTRUE(all.equal(baryRatInitRef(fEC, -pi, pi, 3L, 3L),
                              chebNodes2(8L, -pi, pi))))
expect_error(minimaxApprox(fEC, -pi, pi, c(3L, 3L), basis = "b"),
             "degenerate or defective")
expect_error(minimaxApprox(fEC, -pi, pi, c(14L, 13L), basis = "b"),
             "degenerate or defective")

# E5 re-baseline: pre-E5 the (8, 8) exchange hit a mid-iteration
# no-pole-free stop in this container (the suite's in-loop dispatch
# coverage; M5P2's platform-fragility flag). The E5 exchange converges it
# outright -- E = 2.897293e-7 with dense-grid ratio 1.0000, a global
# certificate, warning-free. Accept either honest outcome; the in-loop
# dispatch line's coverage disposition moves to the NOCOV registry per the
# documented M5P2 fallback.
r88 <- tryCatch(suppressWarnings(minimaxApprox(fEC, -pi, pi, c(8L, 8L),
                                               basis = "b")),
                error = function(e) e)
if (inherits(r88, "error")) {
  expect_true(grepl("degenerate or defective|has a zero at",
                    conditionMessage(r88)))
} else {
  g88 <- seq(-pi, pi, length.out = 2e5L)
  expect_true(max(abs(minimaxEval(g88, r88) - fEC(g88))) / r88$ExpErr <=
                minimaxApprox:::REFLOCALTOL || isTRUE(r88$Warning))
}

# x^4 (2,2): even-monomial NON-NORMAL family. The true minimax equioscillates
# on SEVEN points (the reduced (1,1)-in-x^2 problem's four alternations pull
# back to seven by evenness), which a fixed m+n+2 = 6 reference cannot hold,
# so the exchange has multiple fixed points and the outcome is
# platform-(BLAS-)bistable -- measured: reviewer container detects an
# interior denominator zero and stops; HOMEDESKTOP converges to a
# sub-optimal 6-point-leveled fixed point (reported E = 0.02554 vs true
# E = E_{1,1}(x^2 on [0,1]) = 0.0295085; that result's true sup is 0.0381),
# now flagged by the reference-local certificate warning. Accept either
# honest outcome. NOTE: which branch runs -- and therefore which lines the
# refLocal warning coverage rides on -- depends on the local basin.
r4 <- tryCatch(sW(minimaxApprox(function(x) x ^ 4, -1, 1, c(2L, 2L),
                                basis = "b")),
               error = function(e) e)
if (inherits(r4, "error")) {
  expect_true(grepl("has a zero at|degenerate or defective",
                    conditionMessage(r4)))
} else {
  # The reference alternant is a de la Vallee Poussin LOWER bound on the
  # true minimax error, so a converged ExpErr can never exceed it.
  expect_true(r4$ExpErr <= 0.0295086)

  # E5 re-baseline: a platform where the redesigned exchange converges this
  # non-normal case globally returns it warning-free WITH its dense-grid
  # certificate (Warning FALSE at this magnitude implies the certificate
  # passed, which forces ExpErr at the true minimax 0.0295085);
  # reference-local convergence must still warn.
  g4 <- seq(-1, 1, length.out = 2e5L)
  r4Ratio <- max(abs(minimaxEval(g4, r4) - g4 ^ 4)) / r4$ExpErr
  expect_true(isTRUE(r4$Warning) || r4Ratio <= minimaxApprox:::REFLOCALTOL)
}

# exp(cos(x)) (4,4): same non-normal family, bistable the OTHER way --
# measured: reviewer container converges at E = 1.560e-3 with grid-sup/E =
# 1.11 and the reference-local warning; HOMEDESKTOP stops with the
# no-pole-free-eigenpair error. Accept either honest outcome.
rec <- tryCatch(sW(minimaxApprox(fEC, -pi, pi, c(4L, 4L), basis = "b")),
                error = function(e) e)
if (inherits(rec, "error")) {
  expect_true(grepl("degenerate or defective|has a zero at",
                    conditionMessage(rec)))
} else {
  # E5 re-baseline: the E5 exchange converges this case globally --
  # E = 1.566953e-3 (>= the pre-E5 reference-local lower bound 1.560e-3, as
  # de la Vallee Poussin requires) with dense-grid ratio 1.0000 and no
  # warning. A converged result is acceptable silent ONLY with its global
  # certificate; reference-local convergence must still warn.
  gec <- seq(-pi, pi, length.out = 2e5L)
  recRatio <- max(abs(minimaxEval(gec, rec) - fEC(gec))) / rec$ExpErr
  expect_true(isTRUE(rec$Warning) ||
                recRatio <= minimaxApprox:::REFLOCALTOL)
}

# exp(cos(x)) (2,2): platform-STABLE pre-loop pole stop -- the deterministic
# AAA-fallback initial reference and first solve place a denominator zero at
# -0.319964 on both the reviewer container and HOMEDESKTOP (verified
# identical to six digits), so this expectation carries the pole-stop
# coverage on every platform.
expect_error(minimaxApprox(fEC, -pi, pi, c(2L, 2L), basis = "b"),
             "has a zero at")

# exp(cos(x)) (5,4): the POSITIVE side of the non-normal family, and the
# workaround the reference-local warning recommends. For an even function
# the (4,4)-block best is even and its 6-point alternant in t = x^2 pulls
# back to ELEVEN x-extrema; a (5,4) request supplies exactly m + n + 2 = 11
# reference slots, so the full alternant fits and the exchange certifies
# globally. Triple oracle for E = 1.5669528e-3: barycentric on the reviewer
# container (1.566952757e-3), barycentric on HOMEDESKTOP (1.566953e-3), and
# the CLASSICAL path on the parity-reduced problem exp(cos(sqrt(t))) type
# (2,2) on [0, pi^2] (1.566952756e-3). The recovered rational must also be
# (numerically) even.
r54 <- minimaxApprox(fEC, -pi, pi, c(5L, 4L), basis = "b")
expect_equal(r54$ExpErr, 1.5669528e-3, tolerance = 1e-5)
expect_false(r54$Warning)
expect_true(max(abs(c(r54$a[c(2L, 4L, 6L)], r54$b[c(2L, 4L)]))) < 1e-9)

# fn identically zero: the norm fallback (normf <- 1) must engage and the
# request end in the documented degenerate error, not a low-level failure.
expect_error(minimaxApprox(function(x) x * 0, -1, 1, c(1L, 1L), basis = "b"),
             "degenerate or defective")

# Opts-driven exits. maxiter break: miniter > maxiter blocks isConverged and
# the default conviter (30) blocks the stall exit, so the loop must run to
# maxiter exactly. Stall exit: miniter blocks isConverged while conviter = 2
# lets isUnchanging stop the loop well before maxiter.
mx <- sW(minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b",
                       opts = list(maxiter = 5L, miniter = 6L)))
expect_identical(mx$iterations, 5L)
expect_true(mx$Warning)
st <- sW(minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b",
                       opts = list(miniter = 50L, maxiter = 60L,
                                   conviter = 2L)))
expect_true(st$iterations < 40L)
expect_true(st$Warning)

# showProgress prints the per-iteration line.
expect_message(sW(minimaxApprox(exp, -1, 1, c(2L, 2L), basis = "b",
                                opts = list(maxiter = 3L, miniter = 4L,
                                            showProgress = TRUE))),
               "i: 1 E: ")
expect_error(baryRatFail("nan", 3L, 2L), "error curve is undefined")
