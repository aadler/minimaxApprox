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

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

/////////////////////////// Compensated Horner /////////////////////////////////

// Based on Langlois et al. (2006)
// https://drops.dagstuhl.de/opus/volltexte/2006/442/

// This is the y component of twoSum; the x component is the sum itself. These
// variables need to be initialized as "volatile" as we specifically DO NOT WANT
// the compiler to optimize them out (z to a for example). The entire point of
// EFT algorithms is to capture the floating-point error as best possible!
double twoSumy(double a, double b) {
  volatile double x = a + b;
  volatile double z = x - a;
  return((a - (x - z)) + (b - z));
}

// This is the y component of twoProdFMA; the x component is the product itself.
double twoProdFMAy(double a, double b) {
  double x = a * b;
  return(fma(a, b, -x));
}

// Compensated Horner method combining the Error-Free Transformation component
// and the HornerSum component.
extern SEXP compHorner_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);
  const int nm1 = n - 1;   // Used often

  double *px = REAL(x);
  double *pa = REAL(a);

  // ret will be the final value returned to R. In the EFTHorner portion of the
  // algorithm it replaces s0. In the HornerSum portion it replaces r0.
  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *pret = REAL(ret);

  // n must be at least 1 since degree must be >= 0. Therefore initialize ret
  // with the last value of a. If n == 1 then there is no need to calculate pi,
  // sig, or correction.
  for (int i = 0; i < m; ++i) {
    pret[i] = pa[nm1];
    }

  // If n > 1, we need to use the Compensated Horner scheme to evaluate the
  // polynomial. The algorithm will follow Langlois et al.(2006), which as a
  // Horner method, starts at the end and works backwards; thus j.
  //
  // By traversing i in the outer loop, each "x" is complete after each inner
  // loop. The inner loop calculates the "standard" return and the correction
  // using pi and sigma. Since the pi and sigma are unique to the i/j
  // combination, both EFT and Horner sum can be combined in the inner loop.
  // Once the inner loop finishes, so too is the correction for that x, so it
  // applied at the end of the outer loop. Now only one nested loop is needed.

  if (n > 1) {
    double Ax;
    double pi;
    double sig;
    double correction;

    // Error-Free-Transformation (EFT) Horner AND Horner Sum portion of Langlois
    // et al. (2006).
    for (int i = 0; i < m; ++i) {
      correction = 0.0;
      for (int j = nm1; j-- > 0; ) {
        // EFT
        Ax = pret[i] * px[i];
        pi = twoProdFMAy(pret[i], px[i]);
        pret[i] = Ax + pa[j];
        sig = twoSumy(Ax, pa[j]);
        // Horner Sum correction
        correction *= px[i];
        correction += pi + sig;
      }
      pret[i] += correction;
    }
  }

  UNPROTECT(1);
  return(ret);
}

/////////////////////////// Chebyshev Functions/////////////////////////////////

// chebPoly never used directly, so removed and folded into chebMat and chebCalc

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  // Traverse across columns and calculate coefficient for fixed power 0 < j < k
  // for each x.
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
// forth, it is not good coding practice to have two seperate functions doing
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

static const R_CallMethodDef CallEntries[] = {
  {"compHorner_c",    (DL_FUNC) &compHorner_c,  2},
  {"chebMat_c",       (DL_FUNC) &chebMat_c,     2},
  {"chebCalc_c",      (DL_FUNC) &chebCalc_c,    2},
  {NULL,              NULL,                     0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
