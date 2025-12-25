// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>

#include "compHorner.h"

// Based on Langlois et al. (2006)
// https://drops.dagstuhl.de/opus/volltexte/2006/442/

// This is the y component of twoSum; the x component is the sum itself. These
// variables need to be initialized as "volatile" as we specifically DO NOT WANT
// the compiler to optimize them out (z to a for example). The entire point of
// EFT algorithms is to capture the floating-point error as best possible!
volatile long double twoSumy(long double a, long double b) {
  volatile long double x = a + b;
  volatile long double z = x - a;
  return((a - (x - z)) + (b - z));
}

// This is the y component of twoProdFMA; the x component is the product itself.
volatile long double twoProdFMAy(long double a, long double b) {
  long double x = a * b;
  return(fmal(a, b, -x));
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

  // Create auxiliary long double variables (ending in l) to leverage internal
  // functions at higher precisions and only demote to double when populating
  // the SEXP return object. (AA: 2025-12-25)

  if (n > 1) {
    long double Ax;
    long double pi;
    long double sigma;
    long double correction;
    long double pretl;
    long double pxl;
    long double pal;

    for (int i = 0; i < m; ++i) {
      correction = 0.0L;
      pretl = (long double)pret[i];
      pxl = (long double)px[i];
      for (int j = nm1; j-- > 0; ) {
        // Error-Free-Transformation (EFT) Horner
        pal = (long double)pa[j];
        Ax = pretl * pxl;
        pi = twoProdFMAy(pretl, pxl);
        pretl = Ax + pal;
        sigma = twoSumy(Ax, pal);
        // Horner Sum correction
        correction *= pxl;
        correction += pi + sigma;
      }
      pretl += correction;
      pret[i] = (double)pretl;
    }
  }

  UNPROTECT(1);
  return(ret);
}
