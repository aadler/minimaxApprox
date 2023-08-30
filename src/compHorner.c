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
  volatile double x = a + b;
  volatile double z = x - a;
  return((a - (x - z)) + (b - z));
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
  const int nm1 = n - 1;   // Used often

  // Per Hadley's suggestion to create helper pointer variable once
  // http://adv-r.had.co.nz/C-interface.html#c-vectors - Accessing vector data
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
  // Horner method, starts at the end and works backwards which explains the j
  // decrementation.
  //
  // Also, looking at the algorithm, one doesn't need the entire matrix---just
  // the "last" value. And one can simply overwrite old with new since each
  // cell, once calculated, is never called on for a new cell. Since in C we are
  // working one cell at a time anyway, we can read the old and write the new to
  // it directly.
  //
  // By reversing the order of the loop variables and traversing the i's first,
  // each "x" is complete after an outer loop. The inner loop calculates the
  // "standard" return and the correction using pi and sigma. Since the pi and
  // sigma are unique to the i/j combination, both EFT and Horner summ can be
  // combined in inner loop. Once the inner loop finishes, so is the correction
  // for that x, so it applied at the end of the outer loop. Now only one nested
  // loop is needed.
  //
  // (AA: 2023-08-29)

  if (n > 1) {
    // x element of twoProdFMA used more than once so define as variable
    volatile double Ax;
    volatile double pi;
    volatile double sig;
    double correction[m];
    memset(correction, 0, m * sizeof(double));

    // Error-Free-Transformation (EFT) Horner AND Horner Sum portion of Langlois
    // et al. (2006).
    for (int i = 0; i < m; ++i) {
      for (int j = nm1; j-- > 0; ) {
        // EFT
        Ax = pret[i] * px[i];
        pi = twoPrody(pret[i], px[i]);
        pret[i] = Ax + pa[j];
        sig = twoSumy(Ax, pa[j]);
        // Horner Sum
        correction[i] *= px[i];
        correction[i] += pi + sig;
      }
      // Now apply correction in outer loop
      pret[i] += correction[i];
    }
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
