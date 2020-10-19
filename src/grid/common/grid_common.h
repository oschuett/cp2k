/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_COMMON_H
#define GRID_COMMON_H

/*******************************************************************************
 * \brief Number of Cartesian orbitals up to given angular momentum quantum.
 * \author Ole Schuett
 ******************************************************************************/
static const int ncoset[] = {1,  // l=0
                             4,  // l=1
                             10, // l=2 ...
                             20,  35,  56,  84,  120, 165, 220,  286,
                             364, 455, 560, 680, 816, 969, 1140, 1330};

/*******************************************************************************
 * \brief Tabulation of the Factorial function, e.g. 5! = fac[5] = 120.
 * \author Ole Schuett
 ******************************************************************************/
static const double fac[] = {
    0.10000000000000000000E+01, 0.10000000000000000000E+01,
    0.20000000000000000000E+01, 0.60000000000000000000E+01,
    0.24000000000000000000E+02, 0.12000000000000000000E+03,
    0.72000000000000000000E+03, 0.50400000000000000000E+04,
    0.40320000000000000000E+05, 0.36288000000000000000E+06,
    0.36288000000000000000E+07, 0.39916800000000000000E+08,
    0.47900160000000000000E+09, 0.62270208000000000000E+10,
    0.87178291200000000000E+11, 0.13076743680000000000E+13,
    0.20922789888000000000E+14, 0.35568742809600000000E+15,
    0.64023737057280000000E+16, 0.12164510040883200000E+18,
    0.24329020081766400000E+19, 0.51090942171709440000E+20,
    0.11240007277776076800E+22, 0.25852016738884976640E+23,
    0.62044840173323943936E+24, 0.15511210043330985984E+26,
    0.40329146112660563558E+27, 0.10888869450418352161E+29,
    0.30488834461171386050E+30, 0.88417619937397019545E+31,
    0.26525285981219105864E+33};

/*******************************************************************************
 * \brief Maps three angular momentum components to a single zero based index.
 * \author Ole Schuett
 ******************************************************************************/
static inline int coset(int lx, int ly, int lz) {
  const int l = lx + ly + lz;
  if (l == 0) {
    return 0;
  } else {
    return ncoset[l - 1] + ((l - lx) * (l - lx + 1)) / 2 + lz;
  }
}

/*******************************************************************************
 * \brief Minimum and maximum for integers (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imin(int x, int y) { return (x < y ? x : y); }
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Equivalent of Fortran's MODULO, which always return a positive number.
 *        https://gcc.gnu.org/onlinedocs/gfortran/MODULO.html
 * \author Ole Schuett
 ******************************************************************************/
static inline int modulo(int a, int m) { return ((a % m + m) % m); }

/*******************************************************************************
 * \brief Prototype for BLAS DGEMM subroutine.
 * \author Ole Schuett
 ******************************************************************************/
void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta, double *c,
            const int *ldc);

/*******************************************************************************
 * \brief Convenient wrapper for BLAS DGEMM subroutine.
 * \author Ole Schuett
 ******************************************************************************/
static inline void blas_dgemm(const char transa, const char transb, const int m,
                              const int n, const int k, const double alpha,
                              const double *a, const int lda, const double *b,
                              const int ldb, const double beta, double *c,
                              const int ldc) {
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
         &ldc);
}

/*******************************************************************************
 * \brief Prototype for BLAS DGER subroutine.
 * \author Ole Schuett
 ******************************************************************************/
void dger_(const int *m, const int *n, const double *alpha, const double *x,
           const int *incx, const double *y, const int *incy, const double *a,
           const int *lda);

/*******************************************************************************
 * \brief Convenient wrapper for BLAS DGER subroutine.
 * \author Ole Schuett
 ******************************************************************************/
static inline void blas_dger(const int m, const int n, const double alpha,
                             const double *x, const int incx, const double *y,
                             const int incy, const double *a, const int lda) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

/*******************************************************************************
 * \brief Prototype for BLAS DGEMV subroutine.
 * \author Ole Schuett
 ******************************************************************************/
void dgemv_(const char *transa, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, const double *x, const int *incx,
            const double *beta, const double *y, const int *incy);

/*******************************************************************************
 * \brief Convenient wrapper for BLAS DGEMV subroutine.
 * \author Ole Schuett
 ******************************************************************************/
static inline void blas_dgemv(const char transa, const int m, const int n,
                              const double alpha, const double *a,
                              const int lda, const double *x, const int incx,
                              const double beta, const double *y,
                              const int incy) {
  dgemv_(&transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

#endif // GRID_COMMON_H

// EOF
