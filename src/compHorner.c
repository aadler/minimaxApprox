// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// Most of the _c functions were used while developing the C code. In actuality,
// they should never be called from R as it is faster to use C as much as
// possible. Since there is no good way to pass list results to C functions, the
// core functionality has been rewritten as non-SEXP C and only the "outermost"
// calls will be exposed. The code will be commented/nocoved out, but left in
// the files for pedagogic, debugging, and reminder reasons.
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
  SEXP s0 = PROTECT(allocVector(REALSXP, m));
  SEXP s1 = PROTECT(allocVector(REALSXP, m));
  SEXP piM = PROTECT(allocMatrix(REALSXP, nm1, m));
  SEXP sigM = PROTECT(allocMatrix(REALSXP, nm1, m));
  SEXP val = PROTECT(allocVector(REALSXP, m));
  SEXP ret = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("val"));
  SET_STRING_ELT(names, 1, mkChar("pi"));
  SET_STRING_ELT(names, 2, mkChar("sig"));

  double *px = REAL(x);
  double *pa = REAL(a);
  double *ps0 = REAL(s0);
  double *ps1 = REAL(s1);
  double *pval = REAL(val);
  double *ppiM = REAL(piM);
  double *psigM = REAL(sigM);
  volatile double Ax;

  memset(ps0, 0, m * sizeof(double));
  memset(ps1, 0, m * sizeof(double));
  memset(pval, 0, m * sizeof(double));
  memset(ppiM, 0, m * nm1 * sizeof(double));
  memset(psigM, 0, m * nm1 * sizeof(double));

  if (n > 0) {
    for (int i = 0; i < m; ++i) {
      ps0[i] = pa[nm1];
    }
  }
  if (n > 1) {
  // Key here is to remember R stores matrices as contiguous vectors in
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
        Ax = ps0[i] * px[i];
        ppiM[nm1 * i + j] = twoPrody(ps0[i], px[i]);
        ps1[i] = Ax + pa[j];
        psigM[nm1 * i + j] = twoSumy(Ax, pa[j]);
      }
  // Storing the new as old.
      for (int i = 0; i < m; ++i) {
        ps0[i] = ps1[i];
      }
    }
  }

  SET_VECTOR_ELT(ret, 0, s0);
  SET_VECTOR_ELT(ret, 1, piM);
  SET_VECTOR_ELT(ret, 2, sigM);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(7);
  return(ret);
}

extern SEXP hornerSum_c(SEXP x, SEXP p, SEXP np, SEXP q, SEXP nq) {
  const int m = LENGTH(x);
  int *prows = INTEGER(np);
  int *qrows = INTEGER(nq);

  SEXP ret = PROTECT(allocVector(REALSXP, m));

  double *pret = REAL(ret);
  memset(pret, 0, m * sizeof(double));

  if (prows[0] <= 0) {
    UNPROTECT(1);
    return(ret);
  }

  if (prows[0] != qrows[0]) {
    error("Error polynomials must be of same length.");
  }

  int pcols = LENGTH(p) / prows[0];
  if (pcols != m) {
    error("Polynomials must have same length as x.");
  }

  double *px = REAL(x);
  double *pp = REAL(p);
  double *pq = REAL(q);

  SEXP r0 = PROTECT(allocVector(REALSXP, m));
  SEXP r1 = PROTECT(allocVector(REALSXP, m));
  double *pr0 = REAL(r0);
  double *pr1 = REAL(r1);

  // See above for why this addressing schema is needed.
  for (int i = 0; i < m; ++i) {
    pr0[i] = pp[prows[0] * (i + 1) - 1] +  pq[prows[0] * (i + 1) - 1];
  }

  // See above for why this addressing schema is needed.
  if (prows[0] > 1) {
    int nm1 = prows[0] - 1;
    for (int j = nm1; j-- > 0; ) {
      for (int i = 0; i < m; ++i) {
        pr1[i] = pr0[i] * px[i] + pp[prows[0] * i + j] + pq[prows[0] * i + j];
      }
      for (int i = 0; i < m; ++i) {
        pr0[i] = pr1[i];
      }
    }
  }
  UNPROTECT(3);
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
}
