# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Based on Langlois et al.(2006)
# https://drops.dagstuhl.de/opus/volltexte/2006/442/

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
      A <- twoProd(s[i + 1L, ], x)
      piM[i, ] <- A$y
      B <- twoSum(A$x, a[i])
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

# Vectorized compensatedHorner
compensatedHorner <- function(x, a) {
  EFTH <- eftHorner(x, a)
  cc <- hornerSum(x, EFTH$pi, EFTH$sig)
  EFTH$val + cc
}
