// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <Rmath.h>

#include "minimaxApprox.h"

// Scalar function for Chebyshev polynomial. Makes chebMat a hair slower since
// there may be multiple calls to pow but makes chebCalc faster since it's
// called directly from C. Worth the tradeoff and maintains single-location of
// the Chebyshev calculations for programming safety. This version is faster
// than checking for fabs(x) <= 1 and then a ternary operator to split the else
// branch. (AA: 2024-02-13)

double chebPoly(double x, double k) {
  if (x < -1.0) {
    return(pow(-1.0, k) * cosh(k * acosh(-x)));
  } else if (x <= 1.0) {
    return(cos(k * acos(x)));
  } else {
    return(cosh(k * acosh(x)));
  }
}

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  for (int j = 0; j < n; ++j) {
    int mj = m * j;
    double jj = j;
    for (int i = 0; i < m; ++i) {
      pret[i + mj] = chebPoly(px[i], jj);
    }
  }

  UNPROTECT(1);
  return(ret);
}

extern SEXP chebCalc_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);

  double *px = REAL(x);
  double *pa = REAL(a);

  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *pret = REAL(ret);

  // R's own definitions in DGEMV require the "matrix" to be passed as a
  // pointer to a one-dimensional array. So just like when using SEXP and
  // passing back to R, the matrix is defined as an array and the
  // "row + nrow * col" convention is used.
  double cMat[m * n];

  for (int j = 0; j < n; ++j) {
    int mj = m * j;
    double jj = j;
    for (int i = 0; i < m; ++i) {
      cMat[i + mj] = chebPoly(px[i], jj);
    }
  }

  // Variables needed to conform to R's call of "dgemv".
  char *TR = "N";
  double zero = 0.0;
  double done = 1.0;
  int ONE = 1;
  double (*cM) = cMat;

  F77_CALL(dgemv)(TR, &m, &n, &done, cM, &m, pa, &ONE, &zero, pret, &ONE FCONE);

  UNPROTECT(1);
  return(ret);
}
