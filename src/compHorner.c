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
    double sigma;
    double correction;

    for (int i = 0; i < m; ++i) {
      correction = 0.0;
      for (int j = nm1; j-- > 0; ) {
        // Error-Free-Transformation (EFT) Horner
        Ax = pret[i] * px[i];
        pi = twoProdFMAy(pret[i], px[i]);
        pret[i] = Ax + pa[j];
        sigma = twoSumy(Ax, pa[j]);
        // Horner Sum correction
        correction *= px[i];
        correction += pi + sigma;
      }
      pret[i] += correction;
    }
  }

  UNPROTECT(1);
  return(ret);
}
