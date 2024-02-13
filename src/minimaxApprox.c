// Copyright Avraham Adler (c) 2024
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

// Function to create a matrix of Chebyshev polynomials of order k where k goes
// from 0 to (n-1) for each one of the m entries in the vector x. This  Makes
// chebCalc faster since it's called directly from C instead of Calling chebMat
// inside chebCalc as in earlier development versions. It also maintains the
// single-location of the Chebyshev calculations for programming safety. This
// version is faster than checking for fabs(x) <= 1 and then a ternary operator
// to split the else branch. It also is vectorized internally for a bit of a
// speed up from the prior scalar version. Needs to be passed a pointer to both
// "x" and "ret" which will be the value returned to the functions called from
// R. It is void, because it doesn't return but modifies the ret object in
// place. (AA: 2024-02-13)

void chebPoly(double *ret, double *x, int m, int n) {
  for (int j = 0; j < n; ++j) {
    int mj = m * j;
    double jj = j;
    double monej = pow(-1.0, jj);
    for (int i = 0; i < m; ++i) {
      if (x[i] < -1.0) {
        ret[i + mj] = monej * cosh(jj * acosh(-x[i]));
      } else if (x[i] <= 1.0) {
        ret[i + mj] = cos(jj * acos(x[i]));
      } else {
        ret[i + mj] = cosh(jj * acosh(x[i]));
      }
    }
  }
}

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  chebPoly(pret, px, m, n);

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
  // Need the pointer to cMat for both chebPoly and DGEMV.
  double (*cM) = cMat;

  chebPoly(cM, px, m, n);

  // Remaining variables needed to conform to R's call of "dgemv".
  char *TR = "N";
  double zero = 0.0;
  double done = 1.0;
  int ONE = 1;

  F77_CALL(dgemv)(TR, &m, &n, &done, cM, &m, pa, &ONE, &zero, pret, &ONE FCONE);

  UNPROTECT(1);
  return(ret);
}
