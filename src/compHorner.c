// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// Based on Langlois et al.(2006)
// https://drops.dagstuhl.de/opus/volltexte/2006/442/

// Most of the _c functions were used while developing the C code. In actuality,
// they should never be called from R as it is faster to use C as much as
// possible. Since there is no good way to pass list results to C functions, the
// core functionality has been rewritten as non-SEXP C and only the "outermost"
// calls will be exposed. The code will be commented/nocoved out, but left in
// the files for pedagogic, debugging, and reminder reasons. This is why the
// eftHorner_c and hornerSum_c may have optimzations/efficiencies not found in
// the "component" functions.
//
// Also, many of the variables need to be initialized as "volatile" as we
// specifically DO NOT WANT the compiler to optimize them out. The entire point
// of EFT algorithms is to capture the floating-point error as best possible!
//
// (AA: 2023-08-21)

// # nocov start
extern SEXP horner_c(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);
  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *px = REAL(x);
  double *pa = REAL(a);
  double *pret = REAL(ret);
  for (int i = 0; i < m; ++i) {
    pret[i] = 0;
    for (int j = n; j-- > 0; ) {
      pret[i] = fma(pret[i], px[i], pa[j]);
    }
  }
  UNPROTECT(1);
  return(ret);
}

extern SEXP twoSum_c(SEXP a, SEXP b) {
  const int n = LENGTH(a);
  const int m = LENGTH(b);
  volatile SEXP x = PROTECT(allocVector(REALSXP, n));
  volatile SEXP y = PROTECT(allocVector(REALSXP, n));
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("y"));
  double *pa = REAL(a);
  double *pb = REAL(b);
  double *px = REAL(x);
  double *py = REAL(y);
  volatile double z[n];

  for (int i = 0; i < n; ++i) {
    // In eftHorner, the second term is a singleton. So run this check.
    int j = m == 1 ? 0 : i;
    px[i] = pa[i] + pb[j];
    z[i] = px[i] - pa[i];
    py[i] = (pa[i] - (px[i] - z[i])) + (pb[j] - z[i]);
  }

  SET_VECTOR_ELT(ret, 0, x);
  SET_VECTOR_ELT(ret, 1, y);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(4);
  return(ret);
}

extern SEXP splitA_c(SEXP a) {
// For doubles, q = 53 so r = 27 so magic number 2 ^ r + 1 is 134217729
  const int n = LENGTH(a);
  volatile SEXP x = PROTECT(allocVector(REALSXP, n));
  volatile SEXP y = PROTECT(allocVector(REALSXP, n));
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("h"));
  SET_STRING_ELT(names, 1, mkChar("l"));
  double *pa = REAL(a);
  double *px = REAL(x);
  double *py = REAL(y);
  volatile double z;

  for (int i = 0; i < n; ++i) {
    z = pa[i] * 134217729;
    px[i] = z - (z - pa[i]);
    py[i] = pa[i] - px[i];
  }

  SET_VECTOR_ELT(ret, 0, x);
  SET_VECTOR_ELT(ret, 1, y);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(4);
  return(ret);
}

extern SEXP twoProd_c(SEXP a, SEXP b) {
  const int n = LENGTH(a);
  SEXP x = PROTECT(allocVector(REALSXP, n));
  SEXP y = PROTECT(allocVector(REALSXP, n));
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("y"));
  double *pa = REAL(a);
  double *pb = REAL(b);
  double *px = REAL(x);
  double *py = REAL(y);
  // Below needed if there is no fused-multiply add
  // double Ah, Al, Bh, Bl;

  for (int i = 0; i < n; ++i) {
    px[i] = pa[i] * pb[i];
    py[i] = fma(pa[i], pb[i], -px[i]);
    // Use below instead if there is no fused-multiply add
    // Ah = splithigh(pa[i]);
    // Al = splitlow(pa[i]);
    // Bh = splithigh(pb[i]);
    // Bl = splitlow(pb[i]);
    // py[i] = Al * Bl - (((px[i] - Ah * Bh) - Al * Bh) - Ah * Bl);
  }

  SET_VECTOR_ELT(ret, 0, x);
  SET_VECTOR_ELT(ret, 1, y);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(4);
  return(ret);
}

// Splits are not used if FMA is used.
double splithigh(double a) {
  // For doubles, q = 53 so r = 27 so magic number 2 ^ r + 1 is 134217729
  volatile double r;
  double z = a * 134217729;
  r = z - (z - a);
  return(r);
}

double splitlow(double a) {
  volatile double r;
  r = a - splithigh(a);
  return(r);
}

// # nocov end

double twoSumy(double a, double b) {
  volatile double x, y, z;
  x = a + b;
  z = x - a;
  y = (a - (x - z)) + (b - z);
  return(y);
}

double twoPrody(double a, double b) {
  volatile double x = a * b;
  return(fma(a, b, -x));
}

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
  // The first element of the twoProd is used more than once so save to Ax.
  // See the R code for equivalent.
  volatile double Ax;

  // Since we will return 0 if n = 0---which is technically legal if asking for
  // the best constant estimate or using rational to test polynomial---we must
  // initialize the s0 vector and the pi and sigma matrices with 0.

  memset(ps0, 0, m * sizeof(double));
  memset(ppiM, 0, m * nm1 * sizeof(double));
  memset(psigM, 0, m * nm1 * sizeof(double));

  // If n is at least 1, initialize s0 with the last value of a.
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
  // actually want (numrows * rowIndex + colIndex). See the R code for how it
  // works in R, but it's vectorized, so operations are done for all columns in
  // one row at a time.
  //
  // Also, looking at how twoProd and twoSum work, only need to store A$x
  // equivalent in variable as everything else is a result and not an input too.
  //
  // Lastly, looking at the algorithm, one doesn't need the entire matrix---just
  // the "last" value. When the inner (column) loop finishes, simply store the
  // new results as the old results and then decrement the row counter.
  //
  // (AA: 2023-08-20)
    for (int j = nm1; j-- > 0; ) {
      for (int i = 0; i < m; ++i) {
        // First element of twoProd
        Ax = ps0[i] * px[i];
        // Second element of twoProd
        ppiM[nm1 * i + j] = twoPrody(ps0[i], px[i]);
        // First element of twoSum
        s1[i] = Ax + pa[j];
        // Second element of twoSum
        psigM[nm1 * i + j] = twoSumy(Ax, pa[j]);
      }
  // Storing the new as old.
      for (int i = 0; i < m; ++i) {
        ps0[i] = s1[i];
      }
    }
  }

  SET_VECTOR_ELT(ret, 0, s0);
  SET_VECTOR_ELT(ret, 1, piM);
  SET_VECTOR_ELT(ret, 2, sigM);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(5);
  return(ret);
}

extern SEXP hornerSum_c(SEXP x, SEXP p, SEXP np, SEXP q, SEXP nq) {
  const int m = LENGTH(x);
  const int prows = INTEGER(np)[0];
  const int qrows = INTEGER(nq)[0];
  // Length of the "flattened" matrix. Will be used for error checking.
  const int lp = LENGTH(p);

  // r0 is returned regardles, so declare as an SEXP and initialize here.
  SEXP r0 = PROTECT(allocVector(REALSXP, m));
  double *pr0 = REAL(r0);
  memset(pr0, 0, m * sizeof(double));

  // prows can be 0 if asking for constant. Ignoring size issues in this case, which
  // may be a mistake. Only asking/continuing if prows > 0.
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
    for (int i = 0; i < m; ++i) {
      pr0[i] = pp[prows * (i + 1) - 1] +  pq[prows * (i + 1) - 1];
    }

    // See above for why this addressing schema is needed.
    // Now build "upwards".
    if (prows > 1) {
      int nm1 = prows - 1;
      for (int j = nm1; j-- > 0; ) {
        for (int i = 0; i < m; ++i) {
          r1[i] = pr0[i] * px[i] + (pp[prows * i + j] + pq[prows * i + j]);
        }
        for (int i = 0; i < m; ++i) {
          pr0[i] = r1[i];
        }
      }
    }
  }

  UNPROTECT(1);
  return(r0);
}

static const R_CallMethodDef CallEntries[] = {
  {"horner_c",      (DL_FUNC) &horner_c,    2},
  {"twoSum_c",      (DL_FUNC) &twoSum_c,    2},
  {"splitA_c",      (DL_FUNC) &splitA_c,    1},
  {"twoProd_c",     (DL_FUNC) &twoProd_c,   2},
  {"eftHorner_c",   (DL_FUNC) &eftHorner_c, 2},
  {"hornerSum_c",   (DL_FUNC) &hornerSum_c, 5},
  {NULL,            NULL,                   0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
