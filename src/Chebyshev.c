// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#ifndef  USE_FC_LEN_T
#define  USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif

#include <Rmath.h>

#include "Chebyshev.h"

// Function to create a matrix of Chebyshev polynomials of order k where k goes
// from 0 to (n-1) for each one of the m entries in the vector x. This makes
// chebCalc faster since it's called directly from C instead of calling chebMat
// inside chebCalc as in earlier development versions. It also maintains the
// single-location of the Chebyshev calculations for programming safety. This
// version is faster than checking for fabs(x) <= 1 and then a ternary operator
// to split the else branch. It also is vectorized internally for a bit of a
// speed up from the prior scalar version. Needs to be passed a pointer to both
// "x" and "ret" which will be the value returned to the functions called from
// R. It is void because it doesn't return but modifies the ret object in place.
// (AA: 2024-02-13)
//
// Enough analysis has shown that moving these functions to Fortran does not
// result in enough of a speedup to warrant the added complexity. Stay with C.
// (AA: 2024-05-20)

// F12: m, n, and the mj stride widened from int to R_xlen_t. R_xlen_t is
// ptrdiff_t-width (>= 32 bits, 64-bit on all realistic build targets), so the
// mj = m * j product and the i + mj index can no longer overflow internally.
// The 2^31-cell cap enforced by the two callers below is a deliberate policy
// guard (maintainer decision: hard cap, not long-vector support), not a
// correctness requirement of this function itself.

void chebPolys(double *ret, double *x, R_xlen_t m, R_xlen_t n) {
  // Instead of calling R_pow_di on -1 and j, realize that it's 1 when j = 0 mod
  // 2 and -1 otherwise. As j starts at 0, start with "1" and just keep flipping
  // its sign right before j loops. A bit faster, roughly 1.5%, and elegant.
  int s = 1;
  for (R_xlen_t j = 0; j < n; ++j) {
    R_xlen_t mj = m * j;
    for (R_xlen_t i = 0; i < m; ++i) {
      if (x[i] < -1.0) {
        ret[i + mj] = s * cosh(j * acosh(-x[i]));
      } else if (x[i] <= 1.0) {
        ret[i + mj] = cos(j * acos(x[i]));
      } else {
        ret[i + mj] = cosh(j * acosh(x[i]));
      }
    }
    s *= -1;
  }
}

// F12: shared guard used by both entry points below. m and n arrive as
// R_xlen_t (from XLENGTH / asInteger), so the product is computed in
// R_xlen_t (size_t-width) arithmetic — it cannot itself wrap before the
// comparison runs. The cap is a maintainer policy choice (2^31 cells, i.e.
// no long-vector support for this release; see review doc F12 and the M1
// scope note) rather than a hardware limit, so it is centralized here for
// a single point of future adjustment. Rf_error is used explicitly rather
// than the error() macro: glibc's <error.h> declares an unrelated error()
// with a different signature, and pulling it in by accident silently
// breaks the call.
static void check_cell_cap(R_xlen_t m, R_xlen_t n, const char *fn) {
  if (m > 0 && n > 0 && m > (R_xlen_t) 2147483647 / n) {
    Rf_error("%s: requested matrix (%td x %td) exceeds the supported "
               "2^31-cell limit", fn, (ptrdiff_t) m, (ptrdiff_t) n);
  }
}

extern SEXP chebMat_c(SEXP x, SEXP k) {
  const R_xlen_t m = XLENGTH(x);
  const R_xlen_t n = (R_xlen_t) asInteger(k) + 1;
  check_cell_cap(m, n, "chebMat_c");
  double *px = REAL(x);

  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  double *pret = REAL(ret);

  chebPolys(pret, px, m, n);

  UNPROTECT(1);
  return(ret);
}

extern SEXP chebCalc_c(SEXP x, SEXP a) {
  const R_xlen_t m = XLENGTH(x);
  const R_xlen_t n = XLENGTH(a);
  check_cell_cap(m, n, "chebCalc_c");
  double *px = REAL(x);
  double *pa = REAL(a);

  // F1 fix: the original `double cMat[m * n];` was a C variable-length
  // array allocated on the C call stack. R's default stack (and especially
  // constrained stacks, e.g. some Windows configurations) is only a few MB;
  // for m * n in the low millions this silently overflows the stack and
  // crashes the R process (confirmed: minimaxEval on a 5e6-point grid ->
  // "segfault from C stack overflow"). R_alloc allocates from R's transient
  // heap instead: it survives arbitrarily large m * n (bounded only by the
  // 2^31-cell guard above and available memory), and is automatically
  // reclaimed by R when this call returns or errors -- no matching free(),
  // no leak path, and no risk of the C stack itself blowing up. This also
  // resolves the -Wvla warning (F11's VLA half) since there is no VLA left.
  double *cMat = (double *) R_alloc((size_t) m * (size_t) n, sizeof(double));

  chebPolys(cMat, px, m, n);

  SEXP ret = PROTECT(allocVector(REALSXP, m));
  double *pret = REAL(ret);

  char *TR = "N";
  double d0 = 0.0;
  double d1 = 1.0;
  int i1 = 1;

  // dgemv's m/n/lda arguments are BLAS int*, not R_xlen_t*, so they are
  // narrowed here explicitly. This is safe: the 2^31-cell guard above
  // already ensures m and n individually fit well within INT_MAX whenever
  // the guard passes (their product is capped at 2^31, so each factor is
  // trivially < 2^31 for any n, m >= 1).
  int mi = (int) m;
  int ni = (int) n;

  F77_CALL(dgemv)(TR, &mi, &ni, &d1, cMat, &mi, pa, &i1, &d0, pret, &i1 FCONE);

  UNPROTECT(1);
  return(ret);
}
