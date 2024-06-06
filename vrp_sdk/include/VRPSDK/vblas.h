/**
* Copyright 2021 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
* 
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
* 
*     http://www.apache.org/licenses/LICENSE-2.0
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
**/
/**
 *  @file        vblas.h
 *  @author      Cesar Fuguet-Tortolero
 */
#ifndef __VBLAS_H__
#define __VBLAS_H__

#include <stdint.h>
#include "asm/vpfloat.h"
#include "vmath.h"
#include "common/cache.h"

#ifndef VBLAS_ENABLE_HWPF
#define VBLAS_ENABLE_HWPF 1
#endif
#if VBLAS_ENABLE_HWPF
#include "drivers/hwpf_dcache.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define vblas_likely(x)       __builtin_expect((x),1)
#define vblas_unlikely(x)     __builtin_expect((x),0)

/*
 *  Vector scaling - x = alpha*x
 */
void sscal(int n, float alpha,
           float *x, vpfloat_off_t x_inc, char enable_prefetch);
void dscal(int n, double alpha,
           double *x, vpfloat_off_t x_inc, char enable_prefetch);
void vscal(int precision, int n, const void *alpha, vpfloat_evp_t a_evp,
           void *x, vpfloat_evp_t x_evp, char enable_prefetch);
void cscal (    int n, void * alpha ,
                void *x, vpfloat_off_t x_inc );
void csscal (   int n, float alpha ,
                void *x, vpfloat_off_t x_inc );
void zscal (    int n, void * alpha ,
                void *x, vpfloat_off_t x_inc );
void zsscal (   int n, double alpha ,
                void *x, vpfloat_off_t x_inc );
void vcscal (   int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                void *x, vpfloat_evp_t x_evp );
void vcsscal (  int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                void *x, vpfloat_evp_t x_evp );

/*
 *  Vector copy - y = x
 */
void scopy(int n, const float *x, vpfloat_off_t x_inc,
           float *y, vpfloat_off_t y_inc);
void dcopy(int n, const double *x, vpfloat_off_t x_inc,
           double *y, vpfloat_off_t y_inc);
void vcopy(int n, const void *x, vpfloat_evp_t x_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch);
void vcopy_d_v(int n, const double *x, vpfloat_off_t x_inc,
               void *y, vpfloat_evp_t y_evp, char enable_prefetch);
void vcopy_v_d(int n, const void *x, vpfloat_evp_t x_evp,
               double *y, vpfloat_off_t y_inc, char enable_prefetch);
void ccopy (    int n,
                const void *x, vpfloat_off_t x_inc ,
                void *y, vpfloat_off_t y_inc );
void zcopy (    int n,
                const void *x, vpfloat_off_t x_inc ,
                void *y, vpfloat_off_t y_inc );
void vccopy (   int n,
                const void *x, vpfloat_evp_t x_evp ,
                void *y, vpfloat_evp_t y_evp );

/*
 *  Vector addition with scaling (AXPY) - y = alpha*x + y
 */
void saxpy(int n, float alpha,
           const float *x, vpfloat_off_t x_inc,
           float *y, vpfloat_off_t y_inc, char enable_prefetch);
void daxpy(int n, double alpha,
           const double *x, vpfloat_off_t x_inc,
           double *y, vpfloat_off_t y_inc, char enable_prefetch);
void vaxpy(int precision, int n,
           const void *alpha, vpfloat_evp_t a_evp,
           const void *x, vpfloat_evp_t x_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch);
void caxpy (int n, void * alpha ,
            const void *x, vpfloat_off_t x_inc ,
            void *y, vpfloat_off_t y_inc );
void zaxpy (int n, void * alpha ,
            const void *x, vpfloat_off_t x_inc ,
            void *y, vpfloat_off_t y_inc );
void vcaxpy (   int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                const void *x, vpfloat_evp_t x_evp ,
                void *y, vpfloat_evp_t y_evp );

/*
 *  Scalar vector-vector multiplication (dot product) - x * y
 */
float sdot(int n, const float *x, vpfloat_off_t x_inc,
           const float *y, vpfloat_off_t y_inc, char enable_prefetch);
double ddot(int n, const double *x, vpfloat_off_t x_inc,
            const double *y, vpfloat_off_t y_inc, char enable_prefetch);
void vdot(int precision, int n,
          const void *x, vpfloat_evp_t x_evp,
          const void *y, vpfloat_evp_t y_evp,
          void *r, vpfloat_evp_t r_evp, char enable_prefetch);
void cdotu (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r);
void cdotc (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r);
void zdotu (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r);
void zdotc (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r);
void vcdotu (   int precision , int n,
                const void *x, vpfloat_evp_t x_evp,
                const void *y, vpfloat_evp_t y_evp,
                void *r, vpfloat_evp_t r_evp );
void vcdotc (   int precision , int n,
                const void *x, vpfloat_evp_t x_evp,
                const void *y, vpfloat_evp_t y_evp,
                void *r, vpfloat_evp_t r_evp );

/*
 *  Matrix-vector multiplication
 *
 *  y = (alpha * A * x) + (beta * y)
 *
 *  @param [in]      precision     Working precision of computations
 *  @param [in]      trans         'N' no transposition of A. 'T' transpose A
 *  @param [in]      m             Number of rows in A
 *  @param [in]      n             Number of columns in A
 *  @param [in]      alpha         Pointer to scalar
 *  @param [in]      alpha_evp     Memory environment of alpha scalar
 *  @param [in]      a             Pointer to the A matrix
 *  @param [in]      a_evp         Memory environment of A matrix
 *  @param [in]      lda           Leading dimension of A (allows spacing between rows)
 *  @param [in]      x             Pointer to vector X
 *  @param [in]      x_evp         Memory environment of X vector
 *  @param [in]      b             Pointer to scalar
 *  @param [in]      b_evp         Memory environment of beta scalar
 *  @param [in,out]  y             Pointer to vector Y
 *  @param [in]      y_evp         Memory environment of Y vector
 */
void sgemv(char trans, int m, int n,
           const float alpha,
           const float *a, int lda,
           const float *x, int x_inc,
           const float beta,
           float *y, int y_inc);
void dgemv(char trans, int m, int n,
           const double alpha,
           const double *a, int lda,
           const double *x, int x_inc,
           const double beta,
           double *y, int y_inc);
void dvgemv(int precision, char trans, int m, int n,
            const double alpha,
            const double *a, int lda,
            const void *x, vpfloat_evp_t x_evp,
            const double beta,
            void *y, vpfloat_evp_t y_evp, char enable_prefetch);
void vgemv(int precision, char trans, int m, int n,
           const void *alpha, vpfloat_evp_t alpha_evp,
           const void *a, vpfloat_evp_t a_evp, int lda,
           const void *x, vpfloat_evp_t x_evp,
           const void *beta, vpfloat_evp_t beta_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch);
void cgemv (char trans , int m, int n,
            const void * alpha ,
            const void *a, int lda ,
            const void *x, int x_inc ,
            const void * beta ,
            void *y, int y_inc );
void zgemv (char trans , int m, int n,
            const void * alpha ,
            const void *a, int lda ,
            const void *x, int x_inc ,
            const void * beta ,
            void *y, int y_inc );
void vcgemv (   int precision , char trans , int m, int n,
                const void *alpha , vpfloat_evp_t alpha_evp ,
                const void *a, vpfloat_evp_t a_evp , int lda, 
                const void *x, vpfloat_evp_t x_evp ,
                const void *beta , vpfloat_evp_t beta_evp ,
                void *y, vpfloat_evp_t y_evp );


/**
 *  Triangular system solver
 *
 *  A.X = b
 *
 *  @param [in]      precision     Working precision of computations.
 *  @param [in]      uplo          'L' A with lower triangular shape,
 *                                 'U' A with upper triangular shape.
 *  @param [in]      trans         'N' no transposition of A, 'T' transpose A.
 *  @param [in]      diag          'N' A without unitary diagonal, 'U' A with unitary diagonal.
 *  @param [in]      n             Number of row in A and X.
 *  @param [in]      a             Pointer to the A matrix.
 *  @param [in]      a_evp         Memory environment of A matrix.
 *  @param [in]      lda           Leading dimension of A (allows spacing between rows).
 *  @param [in,out]  x             Pointer to vector X (X contains the vector b in input).
 *  @param [in]      x_evp         Memory environment of X vector.
 */
void vtrsv(int precision, char uplo, char trans, char diag, int n,
           const void * a, vpfloat_evp_t a_evp, int lda,
           void * x, vpfloat_evp_t x_evp);
void dtrsv(char uplo, char trans, char diag, int n,
           const double * a, int lda,
           double * x, int x_inc);
void strsv(char uplo, char trans, char diag, int n,
           const float * a, int lda,
           float * x, int x_inc);

#ifdef __cplusplus
}
#endif

#endif /* __VBLAS_H__ */
