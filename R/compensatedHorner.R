# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Based on Langlois et al.(2006)
# https://drops.dagstuhl.de/opus/volltexte/2006/442/

# Calculation ported to C. Most of the _c functions were used while developing
# the C code. In actuality, they should never be called from R as it is faster
# to use C as much as possible. Since there is no good way to pass list results
# to C functions, the core functionality has been rewritten as non-SEXP C and
# only the "outermost" calls will be exposed. The code will be commented/nocoved
# out, but left in the files for pedagogic, debugging, and reminder reasons.
# (AA: 2023-08-21)

# nocov start
# Vectorized standard horner method (old polyCalc)
horner <- function(x, a) {
  ret <- double(length(x))
  # Using fastest sequence constructor despite it not checking for empty vector
  # as that should not be possible.
  for (i in length(a):1L) { # nolint (see comment above)
    ret <- (ret * x) + a[i]
  }
  ret
}

# Vectorized twoSum
twoSum <- function(a, b) {
  x <- a + b
  z <- x - a
  y <- (a - (x - z)) + (b - z)
  list(x = x, y = y)
}

# Vectorized splitA
splitA <- function(a) {
  # For doubles, q = 53 so r = 27 so magic number 2 ^ r + 1 is 134217729
  z <- a * 134217729
  x <- z - (z - a)
  y <- a - x
  list(h = x, l = y)
}

# Vectorized twoProd
twoProd <- function(a, b) {
  x <- a * b
  A <- splitA(a)
  B <- splitA(b)
  y <- A$l * B$l - (((x - A$h * B$h) - A$l * B$h) - A$h * B$l)
  list(x = x, y = y)
}

# Vectorized eftHorner
eftHorner <- function(x, a) {
  n <- length(a)
  m <- length(x)
  s <- matrix(0, nrow = n, ncol = m)
  piM <- sigM <-  matrix(0, nrow = n - 1L, ncol = m)
  if (n > 0) s[n, ] <- a[n]
  if (n > 1L) {
    for (i in (n - 1L):1) {
      A <- twoProdC(s[i + 1L, ], x)
      piM[i, ] <- A$y
      B <- twoSumC(A$x, a[i])
      s[i, ] <- B$x
      sigM[i, ] <- B$y
    }
  }
  list(val = s[1L, ], pi = piM, sig = sigM)
}

# Vectorized hornerSum
hornerSum <- function(x, p, q) {
  n <- NROW(p)
  if (n <= 0L) {
    return(rep(0, length(x)))
  }
  if (n != NROW(q)) {
    stop("Error polynomials must be of same length.")
  }
  r <- matrix(0, nrow = n, ncol = length(x))
  r[n, ] <- p[n, ] + q[n, ]
  if (n > 1L) {
    for (i in (n - 1L):1) {
      r[i, ] <- r[i + 1L, ] * x + (p[i, ] + q[i, ])
    }
  }
  r[1L, ]
}

# The following C calls are unused in actuality (see .src file) and are also
# kept for debugging etc.
# (AA: 2023-08-20)

twoSumC <- function(a, b) {
  .Call(twoSum_c, as.double(a), as.double(b))
}

splitAC <- function(a) {
  .Call(splitA_c, as.double(a))
}

twoProdC <- function(a, b) {
  .Call(twoProd_c, as.double(a), as.double(b))
}

# nocov end

hornerC <- function(x, a) {
  .Call(horner_c, as.double(x), as.double(a))
}

hornerSumC <- function(x, p, q) {
  .Call(hornerSum_c, as.double(x), as.double(p), NROW(p), as.double(q), NROW(q))
}

eftHornerC <- function(x, a) {
  .Call(eftHorner_c, as.double(x), as.double(a))
}

# Vectorized compensatedHorner
compensatedHorner <- function(x, a) {
  EFTH <- eftHornerC(x, a)
  cc <- hornerSumC(x, EFTH$pi, EFTH$sig)
  EFTH$val + cc
}
