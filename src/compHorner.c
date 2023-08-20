// Copyright Avraham Adler (c) 2023
// SPDX-License-Identifier: MPL-2.0+

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

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
  double z;

  for (int i = 0; i < n; ++i) {
    px[i] = pa[i] + pb[i];
    z = px[i] - pa[i];
    py[i] = (pa[i] - (px[i] - z)) + (pb[i] - z);
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

double splithigh(double a) {
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
  double Ah, Al, Bh, Bl;

  for (int i = 0; i < n; ++i) {
    px[i]= pa[i] * pb[i];
    Ah = splithigh(pa[i]);
    Al = splitlow(pa[i]);
    Bh = splithigh(pb[i]);
    Bl = splitlow(pb[i]);
    py[i] = Al * Bl - (((px[i] - Ah * Bh) - Al * Bh) - Ah * Bl);
  }

  SET_VECTOR_ELT(ret, 0, x);
  SET_VECTOR_ELT(ret, 1, y);
  setAttrib(ret, R_NamesSymbol, names);
  UNPROTECT(4);
  return(ret);
}

static const R_CallMethodDef CallEntries[] = {
  {"horner_c",    (DL_FUNC) &horner_c,  2},
  {"twoSum_c",    (DL_FUNC) &twoSum_c,  2},
  {"splitA_c",    (DL_FUNC) &splitA_c,  1},
  {"twoProd_c",   (DL_FUNC) &twoProd_c, 2},
  {NULL,          NULL,                 0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
