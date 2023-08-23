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
// (AA: 2023-08-22)

// This is the y component of twoSum; the x component is the sum itself.
double twoSumy(double a, double b) {
  volatile double x, y, z;
  x = a + b;
  z = x - a;
  y = (a - (x - z)) + (b - z);
  return(y);
}

// This is the y component of twoProd; the x component is the product itself.
double twoPrody(double a, double b) {
  volatile double x = a * b;
  return(fma(a, b, -x));
}

// Error-Free-Transformation (EFT) Horner method of Langlois et al. (2006)
extern SEXP eftHorner_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);
  const int nm1 = n - 1;
  // s0 is the "older" vector of sums, This is also the "val" that is returned
  // to R so it needs to be an SEXP.
  SEXP s0 = PROTECT(allocVector(REALSXP, m));
  // s1 is the "newer" vector of sums. It does not need to be an SEXP as it
  // never leaves C to go to R.
  double s1[m];
  // This is the pi Matrix which is returned to R.
  SEXP piM = PROTECT(allocMatrix(REALSXP, nm1, m));
  // This is the sigma Matrix which is returned to R.
  SEXP sigM = PROTECT(allocMatrix(REALSXP, nm1, m));
  // This is the return object (an R list).
  SEXP ret = PROTECT(allocVector(VECSXP, 3));
  // This holds the names of the compoments of the return object so it too must
  // be an SEXP.
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  // Set the names.
  SET_STRING_ELT(names, 0, mkChar("val"));
  SET_STRING_ELT(names, 1, mkChar("pi"));
  SET_STRING_ELT(names, 2, mkChar("sig"));

  // Per Hadley's suggestion to create helper pointer variable once
  // http://adv-r.had.co.nz/C-interface.html#c-vectors - Accessing vector data
  double *px = REAL(x);
  double *pa = REAL(a);
  double *ps0 = REAL(s0);
  double *ppiM = REAL(piM);
  double *psigM = REAL(sigM);
  // The first element of twoProd is used more than once so save to Ax.
  volatile double Ax;

  // Since we will return 0 if n = 0---which is technically legal if asking for
  // the best constant estimate or using rational to test polynomial---we must
  // initialize the s0 vector and the pi and sigma matrices with 0.

  memset(ps0, 0, m * sizeof(double));
  memset(ppiM, 0, m * nm1 * sizeof(double));
  memset(psigM, 0, m * nm1 * sizeof(double));

  // If n is at least 1, initialize s0 with the last value of a.
  // Cannot use memset. See https://stackoverflow.com/a/17288891/2726543
  if (n > 0) {
    for (int i = 0; i < m; ++i) {
      ps0[i] = pa[nm1];
    }
  }

  if (n > 1) {
  // If n > 1, we need to use the Compensated Horner scheme to evaluate the
  // polynomial. The algorithm will follow Langlois et al.(2006), which as a
  // Horner method, starts at the end and works backwards.
  //
  // The key here is to remember R stores matrices as contiguous vectors in
  // COLUMN-MAJOR order. So to get to cell [i, j]---assuming C convention that
  // everything starts at 0---you start by iterating over COLUMNS first and then
  // actually want vecindex = (numrows * rowIndex + colIndex). See the R code
  // for how it works in R, but it's vectorized, so operations are done for all
  // columns in one row at a time.
  //
  // Also, looking at how twoProd and twoSum work, only need to store A$x
  // equivalent in variable as everything else is a result and not an input too.
  //
  // Lastly, looking at the algorithm, one doesn't need the entire matrix---just
  // the "last" value. When the inner (column) loop finishes, simply store the
  // new results as the old results and then decrement the row counter.
  //
  // (AA: 2023-08-22)
  int vecidx;
    for (int j = nm1; j-- > 0; ) {
      for (int i = 0; i < m; ++i) {
        // Proper position of value given flattened column-major matrix.
        vecidx = nm1 * i + j;
        // First element of twoProd
        Ax = ps0[i] * px[i];
        // Second element of twoProd
        ppiM[vecidx] = twoPrody(ps0[i], px[i]);
        // First element of twoSum
        s1[i] = Ax + pa[j];
        // Second element of twoSum
        psigM[vecidx] = twoSumy(Ax, pa[j]);
      }
      // Once "row" is finished, copy "new" into "old".
      memcpy(ps0, s1, m * sizeof(double));
    }
  }

  SET_VECTOR_ELT(ret, 0, s0);
  SET_VECTOR_ELT(ret, 1, piM);
  SET_VECTOR_ELT(ret, 2, sigM);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(5);
  return(ret);
}

// Horner sum method of Langlois et al. (2006) to estimate correction.
extern SEXP hornerSum_c(SEXP x, SEXP p, SEXP np, SEXP q) {
  const int m = LENGTH(x);
  const int prows = INTEGER(np)[0];
  // Length of the "flattened" matrix. Will be used for error checking.
  const int lp = LENGTH(p);

  // r0 is returned regardless, so declare as an SEXP and initialize here.
  SEXP r0 = PROTECT(allocVector(REALSXP, m));
  double *pr0 = REAL(r0);
  memset(pr0, 0, m * sizeof(double));

  // prows can be 0 if asking for constant. Ignoring size issues in this case,
  // which may be a mistake. Only asking/continuing if prows > 0.
  if (prows > 0) {
    if (lp != LENGTH(q)) {
      error("Error polynomials must be of same dimension.");
    }

    if (lp / prows != m) {
      error("Polynomials must have same length as x.");
    }

    // If we got here, then prows > 0 and data is of proper size.
    double *px = REAL(x);
    double *pp = REAL(p);
    double *pq = REAL(q);
    // "Newer" sum vector, similar to s1 in eftHorner.
    double r1[m];

    // See above for why this addressing schema is needed.
    // Replace the 0's with the sum of the LAST rows of p and q (pi and Sigma).
    int vecidx;
    for (int i = 0; i < m; ++i) {
      // Proper position of value given flattened column-major matrix.
      vecidx = prows * (i + 1) - 1;
      pr0[i] = pp[vecidx] + pq[vecidx];
    }

    // See above for why this addressing schema is needed.
    // Now build "upwards".
    if (prows > 1) {
      for (int j = prows - 1; j-- > 0; ) {
        for (int i = 0; i < m; ++i) {
          // Proper position of value given flattened column-major matrix.
          vecidx = prows * i + j;
          r1[i] = pr0[i] * px[i] + (pp[vecidx] + pq[vecidx]);
        }
        // Once "row" is finished, copy "new" into "old".
        memcpy(pr0, r1, m * sizeof(double));
      }
    }
  }

  UNPROTECT(1);
  return(r0);
}

static const R_CallMethodDef CallEntries[] = {
  {"eftHorner_c",   (DL_FUNC) &eftHorner_c, 2},
  {"hornerSum_c",   (DL_FUNC) &hornerSum_c, 4},
  {NULL,            NULL,                   0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
