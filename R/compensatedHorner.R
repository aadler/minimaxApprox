# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Based on Langlois et al.(2006)
# https://drops.dagstuhl.de/opus/volltexte/2006/442/

# Vectorized standard Horner method (old polyCalc)
Horner <- function(x, a) {
  ret <- double(length(x))
  # Using fastest sequence constructor despite it not checking for empty vector
  # as that should not be possible.
  for (i in length(a):1L) {
    ret <- (ret * x) + a[i]
  }
  ret
}

# Vectorized TwoSum
TwoSum <- function(a, b) {
  x <- a + b
  z <- x - a
  y <- (a - (x - z)) + (b - z)
  list(x = x, y = y)
}

# Vectorized Split
Split <- function(a) {
  # For doubles, q = 53 so r = 27 so magic number 2 ^ r + 1 is 134217729
  z <- a * 134217729
  x <- z - (z - a)
  y <- a - x
  list(h = x, l = y)
}

# Vectorized TwoProduct
TwoProduct <- function(a, b) {
  x <- a * b
  A <- Split(a)
  B <- Split(b)
  y <- A$l * B$l - (((x - A$h * B$h) - A$l * B$h) - A$h * B$l)
  list(x = x, y = y)
}

# Vectorized EFTHorner
EFTHorner <- function(x, a) {
  n <- length(a)
  m <- length(x)
  s <- matrix(0, nrow = n, ncol = m)
  piM <- sigM <-  matrix(0, nrow = n - 1L, ncol = m)
  if (n > 0) s[n, ] <- a[n]
  if (n > 1L) {
    for (i in (n - 1L):1) {
      A <- TwoProduct(s[i + 1L, ], x)
      piM[i, ] <- A$y
      B <- TwoSum(A$x, a[i])
      s[i, ] <- B$x
      sigM[i, ] <- B$y
    }
  }
  list(val = s[1L, ], pi = piM, sig = sigM)
}

# Vectorized HornerSum
HornerSum <- function(x, p, q) {
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

# Vectorized CompensatedHorner
CompensatedHorner <- function(x, a) {
  EFTH <- EFTHorner(x, a)
  cc <- HornerSum(x, EFTH$pi, EFTH$sig)
  EFTH$val + cc
}
