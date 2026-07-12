# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

# Create the equivalent of a Vandermonde matrix but using Chebyshev polynomials.
# The prior "outer" call on based on the old ChebPoly was ported to C to a
# nested loop. For some reason, calling outer on the C version of chebPoly did
# not work properly. But once the entirety of the matrix build was moved to C,
# and chebPoly removed completely, it became moot.

chebMat <- function(x, k) {
  .Call(chebMat_c, as.double(x), as.double(k))
}

# M6 (F5 Option A). Affine map from the requested range [l, u] to [-1, 1],
# z = 2(x - l)/(u - l) - 1, algebraically rearranged to (2*x - (l+u))/(u-l) to
# keep it a single division. All Chebyshev-basis evaluation now happens on the
# mapped z; the C layer (chebMat_c/chebCalc_c) is untouched and still just
# evaluates T_k at whatever it is given -- this is the single place the map is
# applied, called from every R-level Chebyshev-basis site (evalFunc, polyMat,
# ratMat, interpRescue, checkDenom).
#
# [-1, 1] exactness: when l = -1, u = 1, this reduces to (2x - 0) / 2 = x.
# Both the numerator and the division are by powers of 2, which IEEE 754
# performs exactly (no rounding) for any finite x that does not itself
# overflow/underflow -- so z == x bit-for-bit on the canonical interval. This
# is relied upon, not just expected: it is what makes every existing
# Chebyshev-basis result on [-1, 1] bitwise-identical after this change.
chebMap <- function(x, l, u) {
  (2 * x - (l + u)) / (u - l)
}

# Function to evaluate Chebyshev polynomials and their coefficients. Originally
# was drop(chebMat(x, length(a) - 1L) %*% a). Ported to C and uses R's C
# interface to DGEMV. This is 15%–25% faster than the ported chebMat and %*%.

chebCalc <- function(x, a) {
  .Call(chebCalc_c, as.double(x), as.double(a))
}

# Below based on open-source Netlib function "dconcm" in mathc90.
cheb2mon <- function(a) {
  n <- length(a)
  # The Chebyshev polynomials of order 0 and 1 ARE the monomials x^0 and x^1!
  if (n > 2L) {
    nm2 <- n - 2L
    tp <- 1
    for (j in seq_len(nm2)) {
      # Cannot vectorize next step since already adjusted "a"s affect downstream
      # "a"s as part of the recursion.
      for (i in nm2:j) {
        a[i] <- a[i] - a[i + 2L]
      }
      a[j + 1L] <- a[j + 1L] / 2
      a[j] <- a[j] * tp
      tp <- tp * 2
    }
    a[c(n - 1L, n)] <- a[c(n - 1L, n)] * tp
  }

  a
}

# M6 (F5 Option A). mmA$a/$b are now coefficients of T_k(z), not T_k(x).
# cheb2mon(a) still produces monomial coefficients IN z (unchanged, and
# correctly so -- it is a basis-conversion recursion, agnostic to what the
# free variable is called). aMono/bMono are documented to remain in natural,
# RAW x, so this composes the monomial-in-z polynomial with the affine map
# z = alpha*x + beta (alpha = 2/(u-l), beta = -(u+l)/(u-l)) via the binomial
# theorem: p(alpha*x + beta) = sum_k c_k (alpha*x + beta)^k, expanded and
# re-collected by power of x.
#
# [-1, 1] exactness: an explicit fast path returns c unchanged rather than
# relying on the general formula's floating-point behavior reducing to the
# identity (which it would -- alpha = 1, beta = 0 exactly, so every term with
# j < k contributes beta^(k-j) = 0 -- but the bitwise-identity requirement on
# the canonical interval is strict enough to pin this down explicitly rather
# than trust it).
#
# Conditioning caveat (for documentation): like cheb2mon itself, this is a
# monomial-basis operation and inherits monomial ill-conditioning at high
# degree off [-1, 1] -- the same caveat class already documented for
# cheb2mon/interpRescue's monomial path. aMono/bMono accuracy is bounded by
# this composition, not by the (well-conditioned) mapped-Chebyshev solve that
# produced `a`/`b` in the first place.
composeAffine <- function(c, l, u) {
  if (l == -1 && u == 1) return(c)
  n <- length(c) - 1L
  alpha <- 2 / (u - l)
  beta <- -(u + l) / (u - l)
  d <- numeric(n + 1L)
  for (k in seq_len(n + 1L)) {
    ck <- c[k]
    if (ck == 0) next
    kk <- k - 1L
    for (j in 0:kk) {
      d[j + 1L] <- d[j + 1L] + ck * choose(kk, j) * alpha ^ j * beta ^ (kk - j)
    }
  }
  d
}
