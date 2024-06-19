// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#ifndef  USE_FC_LEN_T
#define  USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif

#include <Rmath.h>

#include "Chebyshev.h"

// Function to create a matrix of Chebyshev polynomials of order k where k goes
// from 0 to (n-1) for each one of the m entries in the vector x. This makes
// chebCalc faster since it's called directly from C instead of calling chebMat
// inside chebCalc as in earlier development versions. It also maintains the
// single-location of the Chebyshev calculations for programming safety. This
// version is faster than checking for fabs(x) <= 1 and then a ternary operator
// to split the else branch. It also is vectorized internally for a bit of a
// speed up from the prior scalar version. Needs to be passed a pointer to both
// "x" and "ret" which will be the value returned to the functions called from
// R. It is void because it doesn't return but modifies the ret object in place.
// (AA: 2024-02-13)
//
// Enough analysis has shown that moving these functions to Fortran does not
// result in enough of a speedup to warrant the added complexity. Stay with C.
// (AA: 2024-05-20)

void chebPolys(double *ret, double *x, int m, int n) {
  // Instead of calling R_pow_di on -1 and j, realize that it's 1 when j = 0 mod
  // 2 and -1 otherwise. As j starts at 0, start with "1" and just keep flipping
  // its sign right before j loops. A bit faster, roughly 1.5%, and elegant.
  int s = 1;
  for (int j = 0; j < n; ++j) {
    int mj = m * j;
    for (int i = 0; i < m; ++i) {
      if (x[i] < -1.0) {
        ret[i + mj] = s * cosh(j * acosh(-x[i]));
      } else if (x[i] <= 1.0) {
        ret[i + mj] = cos(j * acos(x[i]));
      } else {
        ret[i + mj] = cosh(j * acosh(x[i]));
      }
    }
    s *= -1;
  }
}

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  chebPolys(pret, px, m, n);

  UNPROTECT(1);
  return(ret);
}

extern SEXP chebCalc_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);
  double *px = REAL(x);
  double *pa = REAL(a);

  double cMat[m * n];
  double (*pcM) = cMat;

  chebPolys(pcM, px, m, n);

  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *pret = REAL(ret);
  char *TR = "N";
  double d0 = 0.0;
  double d1 = 1.0;
  int i1 = 1;

  F77_CALL(dgemv)(TR, &m, &n, &d1, pcM, &m, pa, &i1, &d0, pret, &i1 FCONE);

  UNPROTECT(1);
  return(ret);
}
