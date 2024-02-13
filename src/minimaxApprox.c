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

// chebPoly never used directly, so removed and folded into chebMat.

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  // Traverse across columns and calculate coefficient for fixed power
  // 0 <= j <= k for each x.
  for (int j = 0; j < n; j++) {
    int mj = m * j;
    double jj = j;
    double monej = pow(-1.0, jj);
    // Chebyshev Polynomial using cos/acos cosh/acosh for each x for fixed "j".
    for (int i = 0; i < m; ++i) {
      if (px[i] < -1.0) {
        pret[i + mj] = monej * cosh(jj * acosh(-px[i]));
      } else if (px[i] < 1.0) {
        pret[i + mj] = cos(jj * acos(px[i]));
      } else {
        pret[i + mj] = cosh(jj * acosh(px[i]));
      }
    }
  }

  UNPROTECT(1);
  return(ret);
}

// While more than 10% more efficient to use all-in-one than to pass back and
// forth, it is not good coding practice to have two separate functions doing
// the same thing, as that leads to bugs when a change is made to one and not
// the other. Therefore, all cheb matricies will be made by chebMat_c and we
// will just live with the overhead of passing them back and forth between C and
// R for chebCalc purposes. Leaving in the commented out code for now for my own
// reference. May remove at any time since not exposed to the API.
// (AA: 2024-02-01)

// extern SEXP chebCalc2_c(SEXP x, SEXP a) {
//   const int m = LENGTH(x);
//   const int n = LENGTH(a);
//
//   double *px = REAL(x);
//   double *pa = REAL(a);
//
//   SEXP ret = PROTECT(allocVector(REALSXP, m));
//   double *pret = REAL(ret);
//
//   // R's own definitions in DGEMV require the "matrix" to be passed as a
//   // pointer to a one-dimensional array. So just like when using SEXP and
//   // passing back to R, the matrix is defined as an array and the
//   // "row + nrow * col" convention is used.
//   double cMat[m * n];
//
//   for (int j = 0; j < n; ++j) {
//     int mj = m * j;
//     double jj = j;
//     double monej = pow(-1.0, jj);
//     for (int i = 0; i < m; ++i) {
//       if (px[i] < -1.0) {
//         cMat[i + mj] = monej * cosh(jj * acosh(-px[i]));
//       } else if (px[i] < 1.0) {
//         cMat[i + mj] = cos(jj * acos(px[i]));
//       } else {
//         cMat[i + mj] = cosh(jj * acosh(px[i]));
//       }
//     }
//   }
//
//   // Variables needed to conform to R's call of "dgemv".
//   char *TR = "N";
//   double zero = 0.0;
//   double done = 1.0;
//   int ONE = 1;
//   double (*cM) = cMat;
//
//   F77_CALL(dgemv)(TR, &m, &n, &done, cM, &m, pa, &ONE, &zero, pret, &ONE FCONE);
//
//   UNPROTECT(1);
//   return(ret);
// }

extern SEXP chebCalc_c(SEXP x, SEXP a) {
  // x is an R matrix built by chebMat_c above.
  const int m = Rf_nrows(x);
  const int n = LENGTH(a);

  double *px = REAL(x);
  double *pa = REAL(a);

  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *pret = REAL(ret);

  // Variables needed to conform to R's call of "dgemv".
  char *TR = "N";
  double zero = 0.0;
  double done = 1.0;
  int ONE = 1;

  F77_CALL(dgemv)(TR, &m, &n, &done, px, &m, pa, &ONE, &zero, pret, &ONE FCONE);

  UNPROTECT(1);
  return(ret);
}
