// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// Based on Langlois et al. (2006)
// https://drops.dagstuhl.de/opus/volltexte/2006/442/
//
// Many of the variables need to be initialized as "volatile" as we specifically
// DO NOT WANT the compiler to optimize them out. The entire point of EFT
// algorithms is to capture the floating-point error as best possible!
//
// (AA: 2023-08-23)

// This is the y component of twoSum; the x component is the sum itself.
double twoSumy(double a, double b) {
  volatile double x, y, z;
  x = a + b;
  z = x - a;
  y = (a - (x - z)) + (b - z);
  return(y);
}

// This is the y component of twoProdFMA; the x component is the product itself.
double twoPrody(double a, double b) {
  volatile double x = a * b;
  return(fma(a, b, -x));
}

// Compensated Horner method combining the Error-Free Transformation component
// and the HornerSum component.
extern SEXP compHorner_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);

  // ret will be the final value returned to R. In the EFTHorner portion of the
  // algorithm it replaces s0. In the HornerSum portion it replaces r0.
  SEXP ret = PROTECT(allocVector(REALSXP, m));

  // Error-Free-Transformation (EFT) Horner portion of Langlois et al. (2006)
  const int nm1 = n - 1;
  double piM[nm1][m];
  double sigM[nm1][m];
  // The first element of twoProd is used more than once so save to Ax.
  volatile double Ax;

  // Per Hadley's suggestion to create helper pointer variable once
  // http://adv-r.had.co.nz/C-interface.html#c-vectors - Accessing vector data
  double *px = REAL(x);
  double *pa = REAL(a);
  double *pret = REAL(ret);

  // Since we will return 0 if n = 0---which is technically legal if asking for
  // the best constant estimate or using rational to test polynomial---we must
  // initialize the ret vector and the pi and sigma matrices with 0.
  memset(pret, 0, m * sizeof(double));
  memset(piM, 0, m * nm1 * sizeof(double));
  memset(sigM, 0, m * nm1 * sizeof(double));

  // If n is at least 1, initialize ret with the last value of a. If n == 1 then
  // there is no need to calculate piM or sigM.
  if (n > 0) {
    for (int i = 0; i < m; ++i) {
      pret[i] = pa[nm1];
    }
  }

  // If n > 1, we need to use the Compensated Horner scheme to evaluate the
  // polynomial. The algorithm will follow Langlois et al.(2006), which as a
  // Horner method, starts at the end and works backwards.
  //
  // Also, looking at the algorithm, one doesn't need the entire matrix---just
  // the "last" value. And one can simply overwrite old with new since each
  // cell, once calculated, is never called on for a new cell. Since in C we are
  // working one cell at a time anyway, we can read the old and write the new to
  // it directly.
  if (n > 1) {
    for (int j = nm1; j-- > 0; ) {
      for (int i = 0; i < m; ++i) {
        Ax = pret[i] * px[i];
        piM[j][i] = twoPrody(pret[i], px[i]);
        pret[i] = Ax + pa[j];
        sigM[j][i] = twoSumy(Ax, pa[j]);
      }
    }
  }

  // Horner Sum portion of Langlois et al. (2006)
  // Use correction vector instead of matrix, as above, and overwrite as we
  // traverse, as above. Initialize and set as 0 which will be the result if a
  // is of length 1 (constant polynomial).
  double correction[m];
  memset(correction, 0, m * sizeof(double));

  // Now build "upwards" if needed.
  if (nm1 > 0) {
    for (int j = nm1; j-- > 0; ) {
      for (int i = 0; i < m; ++i) {
        correction[i] *= px[i];
        correction[i] += piM[j][i] + sigM[j][i];
      }
    }
  }

  // Add correction to value
  for (int i = 0; i < m; ++i) {
    pret[i] += correction[i];
  }

  UNPROTECT(1);
  return(ret);
}

static const R_CallMethodDef CallEntries[] = {
  {"compHorner_c",    (DL_FUNC) &compHorner_c,  2},
  {NULL,              NULL,                     0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
