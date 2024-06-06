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
 *  @file        vblas_vgemv.c
 *  @author      Cesar Fuguet
 */
#include "VRPSDK/vblas.h"
#include "VRPSDK/vutils.h"
#include "VRPSDK/vblas_perfmonitor.h"
#include <stdio.h>

#define VBLAS_PREFETCH_X 1

static void vgemv_rows_chunk_1(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);
static void vgemv_rows_chunk_4(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);
static void vgemv_rows_chunk_8(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);
static void vgemv_cols_chunk_1_trans(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);
static void vgemv_cols_chunk_4_trans(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);
static void vgemv_cols_chunk_8_trans(int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes);

static void dvgemv_rows_chunk_1(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);
//TODO static void dvgemv_rows_chunk_4(
//      int m, int n,
//      const double alpha,
//      const double *a, int lda,
//      const void   *x, int x_bytes,
//      const double beta,
//      void         *y, int y_bytes);
static void dvgemv_rows_chunk_8(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);
static void dvgemv_rows_stride_chunk_8(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);
void dvgemv_cols_chunk_1_trans(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);
void dvgemv_cols_chunk_4_trans(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);
void dvgemv_cols_chunk_8_trans(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void   *x, int x_bytes,
        const double beta,
        void         *y, int y_bytes);


static void vgemv_rows_chunk_1(
        int m, int n,
        const void *alpha,
        const void *a, int a_bytes, int lda,
        const void *x, int x_bytes,
        const void *beta,
        void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr = a_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;

    pcvt_d_p(P24, 0);                       // acc = 0

    //  process the rows one by one
    for (int i = 0; i < m; i++) {
        ple(P0, x_ptr, 0, EVP1);            // x(n)
        x_ptr += x_bytes;
        ple(P1, a_ptr, 0, EVP0);            // A(m)(n)
        a_ptr += a_bytes;

        for (int j = 0; j < (n - 1); j++) {
            pmul(P16, P1, P0, EC0);         // A(m)(n)*x(n)

            ple(P0, x_ptr, 0, EVP1);        // x(n+1)
            ple(P1, a_ptr, 0, EVP0);        // A(m)(n+1)

            x_ptr += x_bytes;
            a_ptr += a_bytes;

            padd(P24, P24, P16, EC0);       // acc += A(m)(n)*x(n)
        }

        //  Last iteration
        //  {{{
        pmul(P16, P1, P0, EC0);             // A(m)(n)*x(n)

        //    Prepare alpha and beta for the next phase
        ple(P0, alpha_ptr, 0, EVP3);        // alpha
        ple(P1, beta_ptr, 0, EVP4);         // beta

        padd(P24, P24, P16, EC0);           // acc += A(m)(n)*x(n)

        a_ptr = a_prev_ptr + a_bytes*lda;
        a_prev_ptr = a_ptr;
        x_ptr = (uintptr_t)x;
        //  }}}

        //  Compute y = beta*y + alpha*acc
        //  {{{
        pmul(P24, P24, P0, EC0);            // acc = alpha*acc
        ple(P2, y_ptr, 0, EVP2);            // y(m) = mem[y(m)]
        pmul(P2, P2, P1, EC0);              // y = beta*y
        padd(P2, P2, P24, EC0);             // y = beta*y + alpha*acc
        pcvt_d_p(P24, 0);                   // acc = 0

        //    Write back resulting element
        pse(P2, y_ptr, 0, EVP2);            // mem(y) = y(m)
        //  }}}

        y_ptr += y_bytes;
    }
}

static void vgemv_cols_chunk_1_trans(int m, int n,
                                     const void *alpha,
                                     const void *a, int a_bytes, int lda,
                                     const void *x, int x_bytes,
                                     const void *beta,
                                     void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;
    const uintptr_t a_row_offset = a_bytes*lda;

    ple(P8, alpha_ptr, 0, EVP3);            // alpha
    ple(P9, beta_ptr, 0, EVP4);             // beta
    for (int j = 0; j < n; j++) {
        pcvt_d_p(P24, 0);                   // acc = 0
        a_prev_ptr = a_ptr;
        for (int i = 0; i < m; i++) {
            ple(P0, a_ptr, 0, EVP0);        // A(n)(m)
            ple(P1, x_ptr, 0, EVP1);        // x(m)
            pmul(P0, P0, P1, EC0);          // A(n)(m)*x(m)
            padd(P24, P24, P0, EC0);        // acc += A(n)(m)*x(m)
            a_ptr += a_row_offset;
            x_ptr += x_bytes;
        }
        pmul(P24, P24, P8, EC0);            // acc = alpha*acc
        a_ptr = a_prev_ptr + a_bytes;
        ple(P25, y_ptr, 0, EVP2);           // y
        pmul(P25, P25, P9, EC0);            // y = beta*y
        padd(P25, P24, P25, EC0);           // y = beta*y + acc
        pse(P25, y_ptr, 0, EVP2);           // mem(y) = y(n)
        x_ptr = (uintptr_t)x;
        y_ptr += y_bytes;
    }
}

static void vgemv_rows_chunk_4(int m, int n,
                               const void *alpha,
                               const void *a, int a_bytes, int lda,
                               const void *x, int x_bytes,
                               const void *beta,
                               void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;
    const uintptr_t a_row_offset = a_bytes*lda;

    pcvt_d_p(P24, 0);                       // acc(0) = 0
    pcvt_d_p(P25, 0);                       // acc(1) = 0
    pcvt_d_p(P26, 0);                       // acc(2) = 0
    pcvt_d_p(P27, 0);                       // acc(3) = 0

    //  process the rows using chunks of 4 rows
    int i;
    for (i = 0; i < m/4; i++) {
        a_prev_ptr = a_ptr;

        ple(P0, x_ptr, 0, EVP1);                  // x(n)
        ple(P1, a_ptr + 0*a_row_offset, 0, EVP0); // A(m+0)(n)
        ple(P2, a_ptr + 1*a_row_offset, 0, EVP0); // A(m+1)(n)
        ple(P3, a_ptr + 2*a_row_offset, 0, EVP0); // A(m+2)(n)
        ple(P4, a_ptr + 3*a_row_offset, 0, EVP0); // A(m+3)(n)

        a_ptr += a_bytes;
        x_ptr += x_bytes;

        for (int j = 0; j < (n - 1); j++) {
            pmul(P5, P1, P0, EC0);                    // A(m+0)(n)*x(n)
            ple(P1, a_ptr + 0*a_row_offset, 0, EVP0); // A(m+0)(n+1)
            pmul(P6, P2, P0, EC0);                    // A(m+1)(n)*x(n)
            ple(P2, a_ptr + 1*a_row_offset, 0, EVP0); // A(m+1)(n+1)
            pmul(P7, P3, P0, EC0);                    // A(m+2)(n)*x(n)
            ple(P3, a_ptr + 2*a_row_offset, 0, EVP0); // A(m+2)(n+1)
            pmul(P8, P4, P0, EC0);                    // A(m+3)(n)*x(n)
            ple(P4, a_ptr + 3*a_row_offset, 0, EVP0); // A(m+3)(n+1)

            ple(P0, x_ptr, 0, EVP1);                  // x(n+1)

            a_ptr += a_bytes;
            x_ptr += x_bytes;

            padd(P24, P24, P5, EC0);                  // acc += A(m+0)(n)*x(n)
            padd(P25, P25, P6, EC0);                  // acc += A(m+1)(n)*x(n)
            padd(P26, P26, P7, EC0);                  // acc += A(m+2)(n)*x(n)
            padd(P27, P27, P8, EC0);                  // acc += A(m+3)(n)*x(n)
        }

        //  Last iteration
        //  {{{
        pmul(P5, P1, P0, EC0);              // A(m+0)(n)*x(n)
        pmul(P6, P2, P0, EC0);              // A(m+1)(n)*x(n)
        pmul(P7, P3, P0, EC0);              // A(m+2)(n)*x(n)
        pmul(P8, P4, P0, EC0);              // A(m+3)(n)*x(n)

        //    Prepare alpha and beta for the next phase
        ple(P0, alpha_ptr, 0, EVP3);        // alpha
        ple(P1, beta_ptr, 0, EVP4);         // beta

        padd(P24, P24, P5, EC0);            // acc += A(m+0)(n)*x(n)
        padd(P25, P25, P6, EC0);            // acc += A(m+1)(n)*x(n)
        padd(P26, P26, P7, EC0);            // acc += A(m+2)(n)*x(n)
        padd(P27, P27, P8, EC0);            // acc += A(m+3)(n)*x(n)

        a_ptr  = a_prev_ptr + 4*a_row_offset;
        x_ptr  = (uintptr_t)x;
        //  }}}

        //  Compute y = beta*y + alpha*acc
        //  {{{
        pmul(P24, P24, P0, EC0);            // acc(0) = acc(0)*alpha
        ple(P2, y_ptr, 0, EVP2);            // y(m+0) = mem[y(m+0)]
        pmul(P25, P25, P0, EC0);            // acc(1) = acc(1)*alpha
        ple(P3, y_ptr, 1, EVP2);            // y(m+1) = mem[y(m+1)]
        pmul(P26, P26, P0, EC0);            // acc(2) = acc(2)*alpha
        ple(P4, y_ptr, 2, EVP2);            // y(m+2) = mem[y(m+2)]
        pmul(P27, P27, P0, EC0);            // acc(3) = acc(3)*alpha
        ple(P5, y_ptr, 3, EVP2);            // y(m+3) = mem[y(m+3)]

        pmul(P2, P2, P1, EC0);              // y(m+0) = y(m+0)*beta
        pmul(P3, P3, P1, EC0);              // y(m+1) = y(m+1)*beta
        pmul(P4, P4, P1, EC0);              // y(m+2) = y(m+2)*beta
        pmul(P5, P5, P1, EC0);              // y(m+3) = y(m+3)*beta

        padd(P2, P2, P24, EC0);             // y(m+0) = beta*y(m+0) + acc(0)
        pcvt_d_p(P24, 0);                   // acc(0) = 0
        padd(P3, P3, P25, EC0);             // y(m+1) = beta*y(m+1) + acc(1)
        pcvt_d_p(P25, 0);                   // acc(1) = 0
        padd(P4, P4, P26, EC0);             // y(m+2) = beta*y(m+2) + acc(2)
        pcvt_d_p(P26, 0);                   // acc(2) = 0
        padd(P5, P5, P27, EC0);             // y(m+3) = beta*y(m+3) + acc(3)
        pcvt_d_p(P27, 0);                   // acc(3) = 0

        //    Write back resulting (partial) vector
        pse(P2, y_ptr, 0, EVP2);            // mem[y(m+0)] = y(m+0)
        pse(P3, y_ptr, 1, EVP2);            // mem[y(m+1)] = y(m+1)
        pse(P4, y_ptr, 2, EVP2);            // mem[y(m+2)] = y(m+2)
        pse(P5, y_ptr, 3, EVP2);            // mem[y(m+3)] = y(m+3)
        //  }}}

        y_ptr += 4*y_bytes;
    }

    //  process the remaining rows (not multiple of 4)
    if ((m % 4) != 0) {
        vgemv_rows_chunk_1(
                m - i*4, n,
                alpha,
                (void*)a_ptr, a_bytes, lda,
                (void*)x_ptr, x_bytes,
                beta,
                (void*)y_ptr, y_bytes);
    }
}

static void vgemv_rows_chunk_8(int m, int n,
                               const void *alpha,
                               const void *a, int a_bytes, int lda,
                               const void *x, int x_bytes,
                               const void *beta,
                               void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;
    const uintptr_t a_row_offset = a_bytes*lda;

    pcvt_d_p(P24, 0);                       // acc(0) = 0
    pcvt_d_p(P25, 0);                       // acc(1) = 0
    pcvt_d_p(P26, 0);                       // acc(2) = 0
    pcvt_d_p(P27, 0);                       // acc(3) = 0
    pcvt_d_p(P28, 0);                       // acc(4) = 0
    pcvt_d_p(P29, 0);                       // acc(5) = 0
    pcvt_d_p(P30, 0);                       // acc(6) = 0
    pcvt_d_p(P31, 0);                       // acc(7) = 0

    //  process the rows using chunks of 8 rows
    int i;
    for (i = 0; i < m/8; i++) {
        a_prev_ptr = a_ptr;

        ple(P0, x_ptr, 0, EVP1);                  // x(n)
        ple(P1, a_ptr + 0*a_row_offset, 0, EVP0); // A(m+0)(n)
        ple(P2, a_ptr + 1*a_row_offset, 0, EVP0); // A(m+1)(n)
        ple(P3, a_ptr + 2*a_row_offset, 0, EVP0); // A(m+2)(n)
        ple(P4, a_ptr + 3*a_row_offset, 0, EVP0); // A(m+3)(n)
        ple(P5, a_ptr + 4*a_row_offset, 0, EVP0); // A(m+4)(n)
        ple(P6, a_ptr + 5*a_row_offset, 0, EVP0); // A(m+5)(n)
        ple(P7, a_ptr + 6*a_row_offset, 0, EVP0); // A(m+6)(n)
        ple(P8, a_ptr + 7*a_row_offset, 0, EVP0); // A(m+7)(n)

        a_ptr += a_bytes;
        x_ptr += x_bytes;

        for (int j = 0; j < (n - 1); j++) {
            pmul(P16, P1, P0, EC0);                   // A(m+0)(n)*x(n)
            ple(P1, a_ptr + 0*a_row_offset, 0, EVP0); // A(m+0)(n+1)
            pmul(P17, P2, P0, EC0);                   // A(m+1)(n)*x(n)
            ple(P2, a_ptr + 1*a_row_offset, 0, EVP0); // A(m+1)(n+1)
            pmul(P18, P3, P0, EC0);                   // A(m+2)(n)*x(n)
            ple(P3, a_ptr + 2*a_row_offset, 0, EVP0); // A(m+2)(n+1)
            pmul(P19, P4, P0, EC0);                   // A(m+3)(n)*x(n)
            ple(P4, a_ptr + 3*a_row_offset, 0, EVP0); // A(m+3)(n+1)
            pmul(P20, P5, P0, EC0);                   // A(m+4)(n)*x(n)
            ple(P5, a_ptr + 4*a_row_offset, 0, EVP0); // A(m+4)(n+1)
            padd(P24, P24, P16, EC0);                 // acc += A(m+0)(n)*x(n)
            padd(P25, P25, P17, EC0);                 // acc += A(m+1)(n)*x(n)
            pmul(P21, P6, P0, EC0);                   // A(m+5)(n)*x(n)
            ple(P6, a_ptr + 5*a_row_offset, 0, EVP0); // A(m+5)(n+1)
            padd(P26, P26, P18, EC0);                 // acc += A(m+2)(n)*x(n)
            padd(P27, P27, P19, EC0);                 // acc += A(m+3)(n)*x(n)
            pmul(P22, P7, P0, EC0);                   // A(m+6)(n)*x(n)
            ple(P7, a_ptr + 6*a_row_offset, 0, EVP0); // A(m+6)(n+1)
            padd(P28, P28, P20, EC0);                 // acc += A(m+4)(n)*x(n)
            padd(P29, P29, P21, EC0);                 // acc += A(m+5)(n)*x(n)
            pmul(P23, P8, P0, EC0);                   // A(m+7)(n)*x(n)
            ple(P8, a_ptr + 7*a_row_offset, 0, EVP0); // A(m+7)(n+1)

            ple(P0, x_ptr, 0, EVP1);                  // x(n+1)
            a_ptr += a_bytes;
            x_ptr += x_bytes;

            padd(P30, P30, P22, EC0);                 // acc += A(m+6)(n)*x(n)
            padd(P31, P31, P23, EC0);                 // acc += A(m+7)(n)*x(n)
        }

        //  Last iteration
        //  {{{
        pmul(P16, P1, P0, EC0);                   // A(m+0)(n)*x(n)
        pmul(P17, P2, P0, EC0);                   // A(m+1)(n)*x(n)
        pmul(P18, P3, P0, EC0);                   // A(m+2)(n)*x(n)
        pmul(P19, P4, P0, EC0);                   // A(m+3)(n)*x(n)
        pmul(P20, P5, P0, EC0);                   // A(m+4)(n)*x(n)
        padd(P24, P24, P16, EC0);                 // acc += A(m+0)(n)*x(n)
        padd(P25, P25, P17, EC0);                 // acc += A(m+1)(n)*x(n)
        pmul(P21, P6, P0, EC0);                   // A(m+5)(n)*x(n)
        padd(P26, P26, P18, EC0);                 // acc += A(m+2)(n)*x(n)
        padd(P27, P27, P19, EC0);                 // acc += A(m+3)(n)*x(n)
        pmul(P22, P7, P0, EC0);                   // A(m+6)(n)*x(n)
        padd(P28, P28, P20, EC0);                 // acc += A(m+4)(n)*x(n)
        padd(P29, P29, P21, EC0);                 // acc += A(m+5)(n)*x(n)
        pmul(P23, P8, P0, EC0);                   // A(m+7)(n)*x(n)

        //    Prepare alpha and beta for the next phase
        ple(P0, alpha_ptr, 0, EVP3);              // alpha
        ple(P1, beta_ptr, 0, EVP4);               // beta

        padd(P30, P30, P22, EC0);                 // acc += A(m+6)(n)*x(n)
        padd(P31, P31, P23, EC0);                 // acc += A(m+7)(n)*x(n)

        a_ptr = a_prev_ptr + 8*a_row_offset;
        x_ptr = (uintptr_t)x;
        //  }}}

        //  Compute y = beta*y + alpha*acc
        //  {{{
        pmul(P24, P24, P0, EC0);            // acc(0) = acc(0)*alpha
        ple(P2, y_ptr, 0, EVP2);            // y(m+0) = mem[y(m+0)]
        pmul(P25, P25, P0, EC0);            // acc(1) = acc(1)*alpha
        ple(P3, y_ptr, 1, EVP2);            // y(m+1) = mem[y(m+1)]
        pmul(P26, P26, P0, EC0);            // acc(2) = acc(2)*alpha
        ple(P4, y_ptr, 2, EVP2);            // y(m+2) = mem[y(m+2)]
        pmul(P27, P27, P0, EC0);            // acc(3) = acc(3)*alpha
        ple(P5, y_ptr, 3, EVP2);            // y(m+3) = mem[y(m+3)]
        pmul(P28, P28, P0, EC0);            // acc(4) = acc(4)*alpha
        ple(P6, y_ptr, 4, EVP2);            // y(m+4) = mem[y(m+4)]
        pmul(P29, P29, P0, EC0);            // acc(5) = acc(5)*alpha
        ple(P7, y_ptr, 5, EVP2);            // y(m+5) = mem[y(m+5)]
        pmul(P30, P30, P0, EC0);            // acc(6) = acc(6)*alpha
        ple(P8, y_ptr, 6, EVP2);            // y(m+6) = mem[y(m+6)]
        pmul(P31, P31, P0, EC0);            // acc(7) = acc(7)*alpha
        ple(P9, y_ptr, 7, EVP2);            // y(m+7) = mem[y(m+7)]

        pmul(P2, P2, P1, EC0);              // y(m+0) = y(m+0)*beta
        pmul(P3, P3, P1, EC0);              // y(m+1) = y(m+1)*beta
        padd(P2, P2, P24, EC0);             // y(m+0) = beta*y(m+0) + acc(0)
        pcvt_d_p(P24, 0);                   // acc(0) = 0
        pmul(P4, P4, P1, EC0);              // y(m+2) = y(m+2)*beta
        padd(P3, P3, P25, EC0);             // y(m+1) = beta*y(m+1) + acc(1)
        pcvt_d_p(P25, 0);                   // acc(1) = 0
        pmul(P5, P5, P1, EC0);              // y(m+3) = y(m+3)*beta
        padd(P4, P4, P26, EC0);             // y(m+2) = beta*y(m+2) + acc(2)
        pcvt_d_p(P26, 0);                   // acc(2) = 0
        pmul(P6, P6, P1, EC0);              // y(m+4) = y(m+4)*beta
        padd(P5, P5, P27, EC0);             // y(m+3) = beta*y(m+3) + acc(3)
        pcvt_d_p(P27, 0);                   // acc(3) = 0
        pmul(P7, P7, P1, EC0);              // y(m+5) = y(m+5)*beta
        padd(P6, P6, P28, EC0);             // y(m+4) = beta*y(m+4) + acc(4)
        pcvt_d_p(P28, 0);                   // acc(4) = 0
        pmul(P8, P8, P1, EC0);              // y(m+6) = y(m+6)*beta
        padd(P7, P7, P29, EC0);             // y(m+5) = beta*y(m+5) + acc(5)
        pcvt_d_p(P29, 0);                   // acc(5) = 0
        pmul(P9, P9, P1, EC0);              // y(m+7) = y(m+7)*beta
        padd(P8, P8, P30, EC0);             // y(m+6) = beta*y(m+6) + acc(6)
        pcvt_d_p(P30, 0);                   // acc(6) = 0
        padd(P9, P9, P31, EC0);             // y(m+7) = beta*y(m+7) + acc(7)
        pcvt_d_p(P31, 0);                   // acc(7) = 0

        //    Write back resulting (partial) vector
        pse(P2, y_ptr, 0, EVP2);            // mem[y(m+0)] = y(m+0)
        pse(P3, y_ptr, 1, EVP2);            // mem[y(m+1)] = y(m+1)
        pse(P4, y_ptr, 2, EVP2);            // mem[y(m+2)] = y(m+2)
        pse(P5, y_ptr, 3, EVP2);            // mem[y(m+3)] = y(m+3)
        pse(P6, y_ptr, 4, EVP2);            // mem[y(m+4)] = y(m+4)
        pse(P7, y_ptr, 5, EVP2);            // mem[y(m+5)] = y(m+5)
        pse(P8, y_ptr, 6, EVP2);            // mem[y(m+6)] = y(m+6)
        pse(P9, y_ptr, 7, EVP2);            // mem[y(m+7)] = y(m+7)
        //  }}}

        y_ptr += 8*y_bytes;
    }

    //  process the remaining rows (not multiple of 8)
    if ((m % 8) != 0) {
        vgemv_rows_chunk_4(m - i*8, n,
                           alpha,
                           (void*)a_ptr, a_bytes, lda,
                           (void*)x_ptr, x_bytes,
                           beta,
                           (void*)y_ptr, y_bytes);
    }
}

static void vgemv_cols_chunk_8_trans(int m, int n,
                                     const void *alpha,
                                     const void *a, int a_bytes, int lda,
                                     const void *x, int x_bytes,
                                     const void *beta,
                                     void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;
    const uintptr_t a_row_offset = a_bytes*lda;

    pcvt_d_p(P24, 0);                       // acc(0) = 0
    pcvt_d_p(P25, 0);                       // acc(1) = 0
    pcvt_d_p(P26, 0);                       // acc(2) = 0
    pcvt_d_p(P27, 0);                       // acc(3) = 0
    pcvt_d_p(P28, 0);                       // acc(4) = 0
    pcvt_d_p(P29, 0);                       // acc(5) = 0
    pcvt_d_p(P30, 0);                       // acc(6) = 0
    pcvt_d_p(P31, 0);                       // acc(7) = 0

    //  process the matrix in chunks of 8 columns
    int j;
    for (j = 0; j < n/8; j++) {
        a_prev_ptr = a_ptr;

        ple(P8, x_ptr, 0, EVP1);                 // x(n)
        x_ptr += x_bytes;

        ple(P0, a_ptr, 0, EVP0);                 // A(n+0)(m)
        ple(P1, a_ptr, 1, EVP0);                 // A(n+1)(m)
        ple(P2, a_ptr, 2, EVP0);                 // A(n+2)(m)
        ple(P3, a_ptr, 3, EVP0);                 // A(n+3)(m)
        ple(P4, a_ptr, 4, EVP0);                 // A(n+4)(m)
        ple(P5, a_ptr, 5, EVP0);                 // A(n+5)(m)
        ple(P6, a_ptr, 6, EVP0);                 // A(n+6)(m)
        ple(P7, a_ptr, 7, EVP0);                 // A(n+7)(m)
        a_ptr += a_row_offset;

        for (int i = 0; i < (m - 1); i++) {
            pmul(P16, P0, P8, EC0);              // A(n+0)(m)*x(n)
            ple(P0, a_ptr, 0, EVP0);             // A(n+0)(m+1)
            pmul(P17, P1, P8, EC0);              // A(n+1)(m)*x(n)
            ple(P1, a_ptr, 1, EVP0);             // A(n+1)(m+1)
            pmul(P18, P2, P8, EC0);              // A(n+2)(m)*x(n)
            ple(P2, a_ptr, 2, EVP0);             // A(n+2)(m+1)
            pmul(P19, P3, P8, EC0);              // A(n+3)(m)*x(n)
            ple(P3, a_ptr, 3, EVP0);             // A(n+3)(m+1)
            padd(P24, P24, P16, EC0);            // acc += A(n+0)(m)*x(n)
            padd(P25, P25, P17, EC0);            // acc += A(n+1)(m)*x(n)
            pmul(P20, P4, P8, EC0);              // A(n+4)(m)*x(n)
            ple(P4, a_ptr, 4, EVP0);             // A(n+4)(m+1)
            padd(P26, P26, P18, EC0);            // acc += A(n+2)(m)*x(n)
            padd(P27, P27, P19, EC0);            // acc += A(n+3)(m)*x(n)
            pmul(P21, P5, P8, EC0);              // A(n+5)(m)*x(n)
            ple(P5, a_ptr, 5, EVP0);             // A(n+5)(m+1)
            padd(P28, P28, P20, EC0);            // acc += A(n+4)(m)*x(n)
            pmul(P22, P6, P8, EC0);              // A(n+6)(m)*x(n)
            ple(P6, a_ptr, 6, EVP0);             // A(n+6)(m+1)
            padd(P29, P29, P21, EC0);            // acc += A(n+5)(m)*x(n)
            pmul(P23, P7, P8, EC0);              // A(n+7)(m)*x(n)
            ple(P8, x_ptr, 0, EVP1);             // x(n+1)
            x_ptr += x_bytes;
            padd(P30, P30, P22, EC0);            // acc += A(n+6)(m)*x(n)
            ple(P7, a_ptr, 7, EVP0);             // A(n+7)(m+1)
            a_ptr += a_row_offset;
            padd(P31, P31, P23, EC0);            // acc += A(n+7)(m)*x(n)
        }

        //  Last iteration
        //  {{{
        pmul(P16, P0, P8, EC0);                  // A(n+0)(m)*x(n)
        ple(P0, y_ptr, 0, EVP2);                 // y(n+0) = mem[y(n+0)]
        a_ptr += a_row_offset;
        pmul(P17, P1, P8, EC0);                  // A(n+1)(m)*x(n)
        ple(P1, y_ptr, 1, EVP2);                 // y(n+1) = mem[y(n+1)]
        x_ptr += x_bytes;
        pmul(P18, P2, P8, EC0);                  // A(n+2)(m)*x(n)
        ple(P2, y_ptr, 2, EVP2);                 // y(n+2) = mem[y(n+2)]
        pmul(P19, P3, P8, EC0);                  // A(n+3)(m)*x(n)
        ple(P3, y_ptr, 3, EVP2);                 // y(n+3) = mem[y(n+3)]
        padd(P24, P24, P16, EC0);                // acc += A(n+0)(m)*x(n)
        padd(P25, P25, P17, EC0);                // acc += A(n+1)(m)*x(n)
        pmul(P20, P4, P8, EC0);                  // A(n+4)(m+0)*x(n)
        ple(P4, y_ptr, 4, EVP2);                 // y(n+0) = mem[y(n+0)]
        padd(P26, P26, P18, EC0);                // acc += A(n+2)(m)*x(n)
        padd(P27, P27, P19, EC0);                // acc += A(n+3)(m)*x(n)
        pmul(P21, P5, P8, EC0);                  // A(n+5)(m+1)*x(n)
        ple(P5, y_ptr, 5, EVP2);                 // y(n+1) = mem[y(n+1)]
        padd(P28, P28, P20, EC0);                // acc += A(n+4)(m)*x(n)
        pmul(P22, P6, P8, EC0);                  // A(n+6)(m+2)*x(n)
        ple(P6, y_ptr, 6, EVP2);                 // y(n+2) = mem[y(n+2)]
        padd(P29, P29, P21, EC0);                // acc += A(n+5)(m)*x(n)
        padd(P30, P30, P22, EC0);                // acc += A(n+6)(m)*x(n)
        pmul(P23, P7, P8, EC0);                  // A(n+7)(m+3)*x(n)
        ple(P8, alpha_ptr, 0, EVP3);             // alpha
        ple(P9, beta_ptr, 0, EVP4);              // beta
        ple(P7, y_ptr, 7, EVP2);                 // y(n+3) = mem[y(n+3)]
        pmul(P24, P24, P8, EC0);                 // acc(0) *= alpha
        padd(P31, P31, P23, EC0);                // acc += A(n+7)(m)*x(n)
        //  }}}

        pmul(P0, P0, P9, EC0);                   // y(n+0) *= beta
        pmul(P25, P25, P8, EC0);                 // acc(1) *= alpha
        pmul(P1, P1, P9, EC0);                   // y(n+1) *= beta
        pmul(P26, P26, P8, EC0);                 // acc(2) *= alpha
        pmul(P2, P2, P9, EC0);                   // y(n+2) *= beta
        padd(P0, P0, P24, EC0);                  // y(n+0) = beta*y(n+0) + alpha*acc(0)
        pcvt_d_p(P24, 0);                        // acc(0) = 0
        pmul(P27, P27, P8, EC0);                 // acc(3) *= alpha
        pmul(P3, P3, P9, EC0);                   // y(n+3) *= beta
        padd(P1, P1, P25, EC0);                  // y(n+1) = beta*y(n+1) + alpha*acc(1)
        pcvt_d_p(P25, 0);                        // acc(1) = 0
        pse(P0, y_ptr, 0, EVP2);                 // mem[y(n+0)] = y(n+0)
        pmul(P28, P28, P8, EC0);                 // acc(4) *= alpha
        pmul(P4, P4, P9, EC0);                   // y(n+4) *= beta
        padd(P2, P2, P26, EC0);                  // y(n+2) = beta*y(n+2) + alpha*acc(2)
        pcvt_d_p(P26, 0);                        // acc(2) = 0
        pse(P1, y_ptr, 1, EVP2);                 // mem[y(n+1)] = y(n+1)
        pmul(P29, P29, P8, EC0);                 // acc(5) *= alpha
        pmul(P5, P5, P9, EC0);                   // y(n+5) *= beta
        padd(P3, P3, P27, EC0);                  // y(n+3) = beta*y(n+3) + alpha*acc(3)
        pcvt_d_p(P27, 0);                        // acc(3) = 0
        pse(P2, y_ptr, 2, EVP2);                 // mem[y(n+2)] = y(n+2)
        pmul(P30, P30, P8, EC0);                 // acc(6) *= alpha
        pmul(P6, P6, P9, EC0);                   // y(n+6) *= beta
        padd(P4, P4, P28, EC0);                  // y(n+4) = beta*y(n+4) + alpha*acc(4)
        pcvt_d_p(P28, 0);                        // acc(4) = 0
        pse(P3, y_ptr, 3, EVP2);                 // mem[y(n+3)] = y(n+3)
        pmul(P31, P31, P8, EC0);                 // acc(7) *= alpha
        pmul(P7, P7, P9, EC0);                   // y(n+7) *= beta
        padd(P5, P5, P29, EC0);                  // y(n+5) = beta*y(n+5) + alpha*acc(5)
        pcvt_d_p(P29, 0);                        // acc(5) = 0
        pse(P4, y_ptr, 4, EVP2);                 // mem[y(n+4)] = y(n+4)
        padd(P6, P6, P30, EC0);                  // y(n+6) = beta*y(n+6) + alpha*acc(6)
        pcvt_d_p(P30, 0);                        // acc(5) = 0
        pse(P5, y_ptr, 5, EVP2);                 // mem[y(n+5)] = y(n+5)
        padd(P7, P7, P31, EC0);                  // y(n+7) = beta*y(n+7) + alpha*acc(7)
        pcvt_d_p(P31, 0);                        // acc(5) = 0
        pse(P6, y_ptr, 6, EVP2);                 // mem[y(n+6)] = y(n+6)
        a_ptr  = a_prev_ptr + 8*a_bytes;
        x_ptr  = (uintptr_t)x;
        pse(P7, y_ptr, 7, EVP2);                 // mem[y(n+7)] = y(n+7)
        y_ptr += 8*y_bytes;
    }

    //  process the remaining columns (not multiple of 8)
    if ((n % 8) != 0) {
        vgemv_cols_chunk_4_trans(m, n - 8*j,
                                 alpha,
                                 (void*)a_ptr, a_bytes, lda,
                                 (void*)x_ptr, x_bytes,
                                 beta,
                                 (void*)y_ptr, y_bytes);
    }
}

static void vgemv_cols_chunk_4_trans(int m, int n,
                                     const void *alpha,
                                     const void *a, int a_bytes, int lda,
                                     const void *x, int x_bytes,
                                     const void *beta,
                                     void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)alpha;
    uintptr_t beta_ptr = (uintptr_t)beta;
    const uintptr_t a_row_offset = a_bytes*lda;

    ple(P8, alpha_ptr, 0, EVP3);            // alpha
    ple(P9, beta_ptr, 0, EVP4);             // beta

    pcvt_d_p(P24, 0);                       // acc(0) = 0
    pcvt_d_p(P25, 0);                       // acc(1) = 0
    pcvt_d_p(P26, 0);                       // acc(2) = 0
    pcvt_d_p(P27, 0);                       // acc(3) = 0

    //  process the columns using chunks of 4 columns
    int j;
    for (j = 0; j < n/4; j++) {
        a_prev_ptr = a_ptr;

        ple(P4, x_ptr, 0, EVP1);                       // x(n)
        ple(P0, a_ptr, 0, EVP0);                       // A(n+0)(m)
        ple(P1, a_ptr, 1, EVP0);                       // A(n+1)(m)
        ple(P2, a_ptr, 2, EVP0);                       // A(n+2)(m)
        ple(P3, a_ptr, 3, EVP0);                       // A(n+3)(m)

        a_ptr += a_row_offset;
        x_ptr += x_bytes;

        for (int i = 0; i < (m - 1); i++) {
            pmul(P16, P0, P4, EC0);                    // A(n+0)(m)*x(n)
            ple(P0, a_ptr, 0, EVP0);                   // A(n+0)(m+1)
            pmul(P17, P1, P4, EC0);                    // A(n+1)(m)*x(n)
            ple(P1, a_ptr, 1, EVP0);                   // A(n+1)(m+1)
            padd(P24, P24, P16, EC0);                  // acc += A(n+0)(m)*x(n)
            pmul(P18, P2, P4, EC0);                    // A(n+2)(m)*x(n)
            ple(P2, a_ptr, 2, EVP0);                   // A(n+2)(m+1)
            padd(P25, P25, P17, EC0);                  // acc += A(n+1)(m)*x(n)
            pmul(P19, P3, P4, EC0);                    // A(n+3)(m)*x(n)
            ple(P4, x_ptr, 0, EVP1);                   // x(n+1)
            ple(P3, a_ptr, 3, EVP0);                   // A(n+3)(m+1)
            padd(P26, P26, P18, EC0);                  // acc += A(n+2)(m)*x(n)
            a_ptr += a_row_offset;
            x_ptr += x_bytes;
            padd(P27, P27, P19, EC0);                  // acc += A(n+3)(m)*x(n)
        }

        pmul(P16, P0, P4, EC0);                        // A(n+0)(m)*x(n)
        ple(P0, y_ptr, 0, EVP2);                       // y(n+0) = mem[y(n+0)]
        pmul(P17, P1, P4, EC0);                        // A(n+1)(m)*x(n)
        ple(P1, y_ptr, 1, EVP2);                       // y(n+1) = mem[y(n+1)]
        pmul(P18, P2, P4, EC0);                        // A(n+2)(m)*x(n)
        ple(P2, y_ptr, 2, EVP2);                       // y(n+2) = mem[y(n+2)]
        padd(P24, P24, P16, EC0);                      // acc(0) += A(n+0)(m)*x(n)
        pmul(P19, P3, P4, EC0);                        // A(n+3)(m)*x(n)
        ple(P3, y_ptr, 3, EVP2);                       // y(n+3) = mem[y(n+3)]
        pmul(P24, P24, P8, EC0);                       // acc(0) *= alpha
        pmul(P0, P0, P9, EC0);                         // y(n+0) = y(n+0)*beta
        padd(P25, P25, P17, EC0);                      // acc(1) += A(n+1)(m)*x(n)
        pmul(P25, P25, P8, EC0);                       // acc(1) *= alpha
        padd(P26, P26, P18, EC0);                      // acc(2) += A(n+2)(m)*x(n)
        pmul(P1, P1, P9, EC0);                         // y(n+1) = y(n+1)*beta
        padd(P27, P27, P19, EC0);                      // acc(3) += A(n+3)(m)*x(n)
        pmul(P26, P26, P8, EC0);                       // acc(2) *= alpha
        pmul(P2, P2, P9, EC0);                         // y(n+2) = y(n+2)*beta
        padd(P0, P0, P24, EC0);                        // y(n+0) = beta*y(n+0) + alpha*acc(0)
        pcvt_d_p(P24, 0);                              // acc(0) = 0
        pmul(P27, P27, P8, EC0);                       // acc(3) *= alpha
        pmul(P3, P3, P9, EC0);                         // y(n+3) = y(n+3)*beta
        padd(P1, P1, P25, EC0);                        // y(n+1) = beta*y(n+1) + alpha*acc(1)
        pse(P0, y_ptr, 0, EVP2);                       // mem[y(n+0)] = y(n+0)
        pcvt_d_p(P25, 0);                              // acc(1) = 0
        padd(P2, P2, P26, EC0);                        // y(n+2) = beta*y(n+2) + alpha*acc(2)
        pse(P1, y_ptr, 1, EVP2);                       // mem[y(n+1)] = y(n+1)
        pcvt_d_p(P26, 0);                              // acc(2) = 0
        padd(P3, P3, P27, EC0);                        // y(n+3) = beta*y(n+3) + alpha*acc(3)
        pse(P2, y_ptr, 2, EVP2);                       // mem[y(n+2)] = y(n+2)
        pcvt_d_p(P27, 0);                              // acc(3) = 0
        a_ptr  = a_prev_ptr + 4*a_bytes;
        x_ptr  = (uintptr_t)x;
        pse(P3, y_ptr, 3, EVP2);                       // mem[y(n+3)] = y(n+3)
        y_ptr += 4*y_bytes;
    }

    //  process the remaining columns (not multiple of 4)
    if ((n % 4) != 0) {
        vgemv_cols_chunk_1_trans(m, n - j*4,
                                 alpha,
                                 (void*)a_ptr, a_bytes, lda,
                                 (void*)x_ptr, x_bytes,
                                 beta,
                                 (void*)y_ptr, y_bytes);
    }
}

void vgemv(int precision, char trans, int m, int n,
           const void *alpha, vpfloat_evp_t alpha_evp,
           const void *a, vpfloat_evp_t a_evp, int lda,
           const void *x, vpfloat_evp_t x_evp,
           const void *beta, vpfloat_evp_t beta_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{

    /*  If constants and matrix are double-precision, call the optimized,
     *  specific version */
    if (vpfloat_evp_is_double(alpha_evp) &&
        vpfloat_evp_is_double(beta_evp) &&
        vpfloat_evp_is_double(a_evp))
    {
        return dvgemv(precision, trans, m, n,
                      *((double*)alpha), (double*)a, lda,
                      x, x_evp,
                      *((double*)beta), y, y_evp, enable_prefetch);
    }

    /* save environment */
    const uint64_t old_ec0 = pger_ec(EC0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_evp2 = pger_evp(EVP2);
    const uint64_t old_evp3 = pger_evp(EVP3);
    const uint64_t old_evp4 = pger_evp(EVP4);

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);
    pser_evp(pack_evp(a_evp.bis, vpfloat_get_rm_mem(),
             a_evp.es, a_evp.stride), EVP0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(),
             x_evp.es, x_evp.stride), EVP1);
    pser_evp(pack_evp(y_evp.bis, vpfloat_get_rm_mem(),
             y_evp.es, y_evp.stride), EVP2);
    pser_evp(pack_evp(alpha_evp.bis, vpfloat_get_rm_mem(),
             alpha_evp.es, alpha_evp.stride), EVP3);
    pser_evp(pack_evp(beta_evp.bis, vpfloat_get_rm_mem(),
             beta_evp.es, beta_evp.stride), EVP4);

    if (trans == 'N') {
        vgemv_rows_chunk_8(
                m, n,
                alpha,
                a, VPFLOAT_SIZEOF(a_evp), lda,
                x, VPFLOAT_SIZEOF(x_evp),
                beta,
                y, VPFLOAT_SIZEOF(y_evp));
    } else {
        vgemv_cols_chunk_8_trans(
                m, n,
                alpha,
                a, VPFLOAT_SIZEOF(a_evp), lda,
                x, VPFLOAT_SIZEOF(x_evp),
                beta,
                y, VPFLOAT_SIZEOF(y_evp));
    }

    /* restore environment */
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_evp(old_evp2, EVP2);
    pser_evp(old_evp3, EVP3);
    pser_evp(old_evp4, EVP4);

    VBLASPERFMONITOR_FUNCTION_END;
}

void sgemv(char trans, int m, int n,
           const float alpha,
           const float *a, int lda,
           const float *x, int x_inc,
           const float beta,
           float *y, int y_inc)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_FLOAT;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_FLOAT;
    y_env.stride = y_inc;

    vgemv(32, trans, m, n,
          (void*)&alpha, VPFLOAT_EVP_FLOAT,
          (void*)a, VPFLOAT_EVP_FLOAT, lda,
          (void*)x, x_env,
          (void*)&beta, VPFLOAT_EVP_FLOAT,
          (void*)y, y_env, 0);

    VBLASPERFMONITOR_FUNCTION_END;
}

void dgemv(char trans, int m, int n,
           const double alpha,
           const double *a, int lda,
           const double *x, int x_inc,
           const double beta,
           double *y, int y_inc)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_DOUBLE;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_DOUBLE;
    y_env.stride = y_inc;

    vgemv(64, trans, m, n,
          (void*)&alpha, VPFLOAT_EVP_DOUBLE,
          (void*)a, VPFLOAT_EVP_DOUBLE, lda,
          (void*)x, x_env,
          (void*)&beta, VPFLOAT_EVP_DOUBLE,
          (void*)y, y_env, 0);

    VBLASPERFMONITOR_FUNCTION_END;

}

static void dvgemv_rows_chunk_8(int m, int n,
                                const double alpha,
                                const double *a, int lda,
                                const void   *x, int x_bytes,
                                const double beta,
                                void         *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    const uintptr_t a_row_offset = sizeof(double)*lda;

#if 0
    static int trace_done = 0;
    if ( trace_done == 0 )  { printf("Prefetch disabled.\n"); trace_done=1; }
#endif

    //  process the rows using chunks of 8 rows
    int i;
    for (i = 0; i < m/8; i++) {
        a_prev_ptr = a_ptr;
        ple(P8, x_ptr, 0, EVP0);                  // x(n)
        x_ptr += x_bytes;
        pld(P0, a_ptr + 0*a_row_offset, 0, EFP0); // A(m+0)(n)
        pcvt_d_p(P24, 0);                         // acc(0) = 0
        pld(P1, a_ptr + 1*a_row_offset, 0, EFP0); // A(m+1)(n)
        pcvt_d_p(P25, 0);                         // acc(1) = 0
        pmul(P16, P0, P8, EC0);                   // A(m+0)(n)*x(n)
        pld(P2, a_ptr + 2*a_row_offset, 0, EFP0); // A(m+2)(n)
        pcvt_d_p(P26, 0);                         // acc(2) = 0
        pld(P3, a_ptr + 3*a_row_offset, 0, EFP0); // A(m+3)(n)
        pcvt_d_p(P27, 0);                         // acc(3) = 0
        pmul(P17, P1, P8, EC0);                   // A(m+1)(n)*x(n)
        pld(P4, a_ptr + 4*a_row_offset, 0, EFP0); // A(m+4)(n)
        pcvt_d_p(P28, 0);                         // acc(4) = 0
        pld(P5, a_ptr + 5*a_row_offset, 0, EFP0); // A(m+5)(n)
        pcvt_d_p(P29, 0);                         // acc(5) = 0
        pmul(P18, P2, P8, EC0);                   // A(m+2)(n)*x(n)
        pld(P6, a_ptr + 6*a_row_offset, 0, EFP0); // A(m+6)(n)
        pcvt_d_p(P30, 0);                         // acc(6) = 0
        pld(P7, a_ptr + 7*a_row_offset, 0, EFP0); // A(m+7)(n)
        pcvt_d_p(P31, 0);                         // acc(7) = 0
        pmul(P19, P3, P8, EC0);                   // A(m+3)(n)*x(n)
        a_ptr += sizeof(double);
        for (int j = 0; j < (n - 1); j++) {
            pld(P0, a_ptr + 0*a_row_offset, 0, EFP0); // A(m+0)(n+1)
            padd(P24, P24, P16, EC0);                 // acc(0) += A(m+0)(n)*x(n)
            pmul(P20, P4, P8, EC0);                   // A(m+4)(n)*x(n)
            pld(P1, a_ptr + 1*a_row_offset, 0, EFP0); // A(m+1)(n+1)
            padd(P25, P25, P17, EC0);                 // acc(1) += A(m+1)(n)*x(n)
            pmul(P21, P5, P8, EC0);                   // A(m+5)(n)*x(n)
            pld(P2, a_ptr + 2*a_row_offset, 0, EFP0); // A(m+2)(n+1)
            padd(P26, P26, P18, EC0);                 // acc(2) += A(m+2)(n)*x(n)
            pmul(P22, P6, P8, EC0);                   // A(m+6)(n)*x(n)
            pld(P3, a_ptr + 3*a_row_offset, 0, EFP0); // A(m+3)(n+1)
            padd(P27, P27, P19, EC0);                 // acc(3) += A(m+3)(n)*x(n)
            pmul(P23, P7, P8, EC0);                   // A(m+7)(n)*x(n)
            ple(P8, x_ptr, 0, EVP0);                  // x(n+1)
            x_ptr += x_bytes;
            pld(P4, a_ptr + 4*a_row_offset, 0, EFP0); // A(m+4)(n+1)
            padd(P28, P28, P20, EC0);                 // acc(4) += A(m+4)(n)*x(n)
            pmul(P16, P0, P8, EC0);                   // A(m+0)(n)*x(n)
            pld(P5, a_ptr + 5*a_row_offset, 0, EFP0); // A(m+5)(n+1)
            padd(P29, P29, P21, EC0);                 // acc(5) += A(m+5)(n)*x(n)
            pmul(P17, P1, P8, EC0);                   // A(m+1)(n)*x(n)
            pld(P6, a_ptr + 6*a_row_offset, 0, EFP0); // A(m+6)(n+1)
            padd(P30, P30, P22, EC0);                 // acc(6) += A(m+6)(n)*x(n)
            pmul(P18, P2, P8, EC0);                   // A(m+2)(n)*x(n)
            pld(P7, a_ptr + 7*a_row_offset, 0, EFP0); // A(m+7)(n+1)
            padd(P31, P31, P23, EC0);                 // acc(7) += A(m+7)(n)*x(n)
            pmul(P19, P3, P8, EC0);                   // A(m+3)(n)*x(n)
            a_ptr += sizeof(double);
        }
        padd(P24, P24, P16, EC0);                     // acc += A(m+0)(n)*x(n)
        ple(P0, y_ptr, 0, EVP1);                      // y(m+0) = mem[y(m+0)]
        pmul(P20, P4, P8, EC0);                       // A(m+4)(n)*x(n)
        padd(P25, P25, P17, EC0);                     // acc += A(m+1)(n)*x(n)
        ple(P1, y_ptr, 1, EVP1);                      // y(m+1) = mem[y(m+1)]
        pmul(P21, P5, P8, EC0);                       // A(m+5)(n)*x(n)
        padd(P26, P26, P18, EC0);                     // acc += A(m+2)(n)*x(n)
        ple(P2, y_ptr, 2, EVP1);                      // y(m+2) = mem[y(m+2)]
        pmul(P22, P6, P8, EC0);                       // A(m+6)(n)*x(n)
        padd(P27, P27, P19, EC0);                     // acc += A(m+3)(n)*x(n)
        ple(P3, y_ptr, 3, EVP1);                      // y(m+3) = mem[y(m+3)]
        pmul(P23, P7, P8, EC0);                       // A(m+7)(n)*x(n)
        padd(P28, P28, P20, EC0);                     // acc += A(m+4)(n)*x(n)
        ple(P4, y_ptr, 4, EVP1);                      // y(m+4) = mem[y(m+4)]
        padd(P29, P29, P21, EC0);                     // acc += A(m+5)(n)*x(n)
        ple(P5, y_ptr, 5, EVP1);                      // y(m+5) = mem[y(m+5)]
        padd(P30, P30, P22, EC0);                     // acc += A(m+6)(n)*x(n)
        ple(P6, y_ptr, 6, EVP1);                      // y(m+6) = mem[y(m+6)]
        padd(P31, P31, P23, EC0);                     // acc += A(m+7)(n)*x(n)
        ple(P7, y_ptr, 7, EVP1);                      // y(m+7) = mem[y(m+7)]
        a_ptr = a_prev_ptr + 8*a_row_offset;
        pld(P9, beta_ptr, 0, EFP0);                   // beta
        x_ptr = (uintptr_t)x;
        pld(P8, alpha_ptr, 0, EFP0);                  // alpha
        pmul(P0, P0, P9, EC0);                        // y(m+0) = y(m+0)*beta
        pmul(P1, P1, P9, EC0);                        // y(m+1) = y(m+1)*beta
        pmul(P2, P2, P9, EC0);                        // y(m+2) = y(m+2)*beta
        pmul(P3, P3, P9, EC0);                        // y(m+3) = y(m+3)*beta
        pmul(P24, P24, P8, EC0);                      // acc(0) = acc(0)*alpha
        pmul(P25, P25, P8, EC0);                      // acc(1) = acc(1)*alpha
        pmul(P26, P26, P8, EC0);                      // acc(2) = acc(2)*alpha
        pmul(P27, P27, P8, EC0);                      // acc(3) = acc(3)*alpha
        padd(P0, P0, P24, EC0);                       // y(m+0) = beta*y(m+0) + acc(0)
        pmul(P28, P28, P8, EC0);                      // acc(4) = acc(4)*alpha
        padd(P1, P1, P25, EC0);                       // y(m+1) = beta*y(m+1) + acc(1)
        pmul(P29, P29, P8, EC0);                      // acc(5) = acc(5)*alpha
        padd(P2, P2, P26, EC0);                       // y(m+2) = beta*y(m+2) + acc(2)
        pmul(P30, P30, P8, EC0);                      // acc(6) = acc(6)*alpha
        padd(P3, P3, P27, EC0);                       // y(m+3) = beta*y(m+3) + acc(3)
        pmul(P31, P31, P8, EC0);                      // acc(7) = acc(7)*alpha
        pse(P0, y_ptr, 0, EVP1);                      // mem[y(m+0)] = y(m+0)
        pse(P1, y_ptr, 1, EVP1);                      // mem[y(m+1)] = y(m+1)
        pse(P2, y_ptr, 2, EVP1);                      // mem[y(m+2)] = y(m+2)
        pse(P3, y_ptr, 3, EVP1);                      // mem[y(m+3)] = y(m+3)
        pmul(P4, P4, P9, EC0);                        // y(m+4) = y(m+4)*beta
        pmul(P5, P5, P9, EC0);                        // y(m+5) = y(m+5)*beta
        pmul(P6, P6, P9, EC0);                        // y(m+6) = y(m+6)*beta
        pmul(P7, P7, P9, EC0);                        // y(m+7) = y(m+7)*beta
        padd(P4, P4, P28, EC0);                       // y(m+4) = beta*y(m+4) + acc(4)
        padd(P5, P5, P29, EC0);                       // y(m+5) = beta*y(m+5) + acc(5)
        padd(P6, P6, P30, EC0);                       // y(m+6) = beta*y(m+6) + acc(6)
        padd(P7, P7, P31, EC0);                       // y(m+7) = beta*y(m+7) + acc(7)
        pse(P4, y_ptr, 4, EVP1);                      // mem[y(m+4)] = y(m+4)
        pse(P5, y_ptr, 5, EVP1);                      // mem[y(m+5)] = y(m+5)
        pse(P6, y_ptr, 6, EVP1);                      // mem[y(m+6)] = y(m+6)
        pse(P7, y_ptr, 7, EVP1);                      // mem[y(m+7)] = y(m+7)
        y_ptr += 8*y_bytes;
    }

    //  process the remaining rows (not multiple of 8)
    if ((m % 8) != 0) {
        dvgemv_rows_chunk_1(m - i*8, n,
                            alpha,
                            (void*)a_ptr, lda,
                            (void*)x_ptr, x_bytes,
                            beta,
                            (void*)y_ptr, y_bytes);
    }
}

static void dvgemv_rows_chunk_1(int m, int n,
                                const double alpha,
                                const double *a, int lda,
                                const void   *x, int x_bytes,
                                const double beta,
                                void         *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr = a_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_offset = a_bytes*lda;

    pcvt_d_p(P24, 0);                       // acc = 0

    //  process the rows one by one
    for (int i = 0; i < m; i++) {
        ple(P0, x_ptr, 0, EVP0);            // x(n)
        x_ptr += x_bytes;
        pld(P1, a_ptr, 0, EFP0);            // A(m)(n)
        a_ptr += a_bytes;

        for (int j = 0; j < (n - 1); j++) {
            pmul(P16, P1, P0, EC0);         // A(m)(n)*x(n)

            ple(P0, x_ptr, 0, EVP0);        // x(n+1)
            pld(P1, a_ptr, 0, EFP0);        // A(m)(n+1)

            x_ptr += x_bytes;
            a_ptr += a_bytes;

            padd(P24, P24, P16, EC0);       // acc += A(m)(n)*x(n)
        }

        //  Last iteration
        //  {{{
        pmul(P16, P1, P0, EC0);             // A(m)(n)*x(n)

        //    Prepare alpha and beta for the next phase
        pld(P0, alpha_ptr, 0, EFP0);        // alpha
        pld(P1, beta_ptr, 0, EFP0);         // beta

        padd(P24, P24, P16, EC0);           // acc += A(m)(n)*x(n)

        a_ptr = a_prev_ptr + a_row_offset;
        a_prev_ptr = a_ptr;
        x_ptr = (uintptr_t)x;
        //  }}}

        //  Compute y = beta*y + alpha*acc
        //  {{{
        pmul(P24, P24, P0, EC0);            // acc = alpha*acc
        ple(P2, y_ptr, 0, EVP1);            // y(m) = mem[y(m)]
        pmul(P2, P2, P1, EC0);              // y = beta*y
        padd(P2, P2, P24, EC0);             // y = beta*y + alpha*acc
        pcvt_d_p(P24, 0);                   // acc = 0

        //    Write back resulting element
        pse(P2, y_ptr, 0, EVP1);            // mem(y) = y(m)
        //  }}}

        y_ptr += y_bytes;
    }
}

static void dvgemv_rows_stride_chunk_8(int m, int n,
                                       const double alpha,
                                       const double *a, int lda,
                                       const void   *x, int x_bytes,
                                       const double beta,
                                       void         *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    const uintptr_t a_row_offset = sizeof(double)*lda;

    //  process the rows using chunks of 8 rows
    int i;
    for (i = 0; i < m/8; i++) {
        a_prev_ptr = a_ptr;
        ple(P8, x_ptr, 0, EVP0);                  // x(n)
        x_ptr += x_bytes;
        pld(P0, a_ptr, 0, EFP0); // A(m+0)(n)
        pcvt_d_p(P24, 0);                         // acc(0) = 0
        pld(P1, a_ptr, 1, EFP0); // A(m+1)(n)
        pcvt_d_p(P25, 0);                         // acc(1) = 0
        pmul(P16, P0, P8, EC0);                   // A(m+0)(n)*x(n)
        pld(P2, a_ptr, 2, EFP0); // A(m+2)(n)
        pcvt_d_p(P26, 0);                         // acc(2) = 0
        pld(P3, a_ptr, 3, EFP0); // A(m+3)(n)
        pcvt_d_p(P27, 0);                         // acc(3) = 0
        pmul(P17, P1, P8, EC0);                   // A(m+1)(n)*x(n)
        pld(P4, a_ptr, 4, EFP0); // A(m+4)(n)
        pcvt_d_p(P28, 0);                         // acc(4) = 0
        pld(P5, a_ptr, 5, EFP0); // A(m+5)(n)
        pcvt_d_p(P29, 0);                         // acc(5) = 0
        pmul(P18, P2, P8, EC0);                   // A(m+2)(n)*x(n)
        pld(P6, a_ptr, 6, EFP0); // A(m+6)(n)
        pcvt_d_p(P30, 0);                         // acc(6) = 0
        pld(P7, a_ptr, 7, EFP0); // A(m+7)(n)
        pcvt_d_p(P31, 0);                         // acc(7) = 0
        pmul(P19, P3, P8, EC0);                   // A(m+3)(n)*x(n)
        a_ptr += sizeof(double);
        for (int j = 0; j < (n - 1); j++) {
            pld(P0, a_ptr, 0, EFP0); // A(m+0)(n+1)
            padd(P24, P24, P16, EC0);                 // acc(0) += A(m+0)(n)*x(n)
            pmul(P20, P4, P8, EC0);                   // A(m+4)(n)*x(n)
            pld(P1, a_ptr, 1, EFP0); // A(m+1)(n+1)
            padd(P25, P25, P17, EC0);                 // acc(1) += A(m+1)(n)*x(n)
            pmul(P21, P5, P8, EC0);                   // A(m+5)(n)*x(n)
            pld(P2, a_ptr, 2, EFP0); // A(m+2)(n+1)
            padd(P26, P26, P18, EC0);                 // acc(2) += A(m+2)(n)*x(n)
            pmul(P22, P6, P8, EC0);                   // A(m+6)(n)*x(n)
            pld(P3, a_ptr, 3, EFP0); // A(m+3)(n+1)
            padd(P27, P27, P19, EC0);                 // acc(3) += A(m+3)(n)*x(n)
            pmul(P23, P7, P8, EC0);                   // A(m+7)(n)*x(n)
            ple(P8, x_ptr, 0, EVP0);                  // x(n+1)
            x_ptr += x_bytes;
            pld(P4, a_ptr, 4, EFP0); // A(m+4)(n+1)
            padd(P28, P28, P20, EC0);                 // acc(4) += A(m+4)(n)*x(n)
            pmul(P16, P0, P8, EC0);                   // A(m+0)(n)*x(n)
            pld(P5, a_ptr, 5, EFP0); // A(m+5)(n+1)
            padd(P29, P29, P21, EC0);                 // acc(5) += A(m+5)(n)*x(n)
            pmul(P17, P1, P8, EC0);                   // A(m+1)(n)*x(n)
            pld(P6, a_ptr, 6, EFP0); // A(m+6)(n+1)
            padd(P30, P30, P22, EC0);                 // acc(6) += A(m+6)(n)*x(n)
            pmul(P18, P2, P8, EC0);                   // A(m+2)(n)*x(n)
            pld(P7, a_ptr, 7, EFP0); // A(m+7)(n+1)
            padd(P31, P31, P23, EC0);                 // acc(7) += A(m+7)(n)*x(n)
            pmul(P19, P3, P8, EC0);                   // A(m+3)(n)*x(n)
            a_ptr += sizeof(double);
        }
        padd(P24, P24, P16, EC0);                     // acc += A(m+0)(n)*x(n)
        ple(P0, y_ptr, 0, EVP1);                      // y(m+0) = mem[y(m+0)]
        pmul(P20, P4, P8, EC0);                       // A(m+4)(n)*x(n)
        padd(P25, P25, P17, EC0);                     // acc += A(m+1)(n)*x(n)
        ple(P1, y_ptr, 1, EVP1);                      // y(m+1) = mem[y(m+1)]
        pmul(P21, P5, P8, EC0);                       // A(m+5)(n)*x(n)
        padd(P26, P26, P18, EC0);                     // acc += A(m+2)(n)*x(n)
        ple(P2, y_ptr, 2, EVP1);                      // y(m+2) = mem[y(m+2)]
        pmul(P22, P6, P8, EC0);                       // A(m+6)(n)*x(n)
        padd(P27, P27, P19, EC0);                     // acc += A(m+3)(n)*x(n)
        ple(P3, y_ptr, 3, EVP1);                      // y(m+3) = mem[y(m+3)]
        pmul(P23, P7, P8, EC0);                       // A(m+7)(n)*x(n)
        padd(P28, P28, P20, EC0);                     // acc += A(m+4)(n)*x(n)
        ple(P4, y_ptr, 4, EVP1);                      // y(m+4) = mem[y(m+4)]
        padd(P29, P29, P21, EC0);                     // acc += A(m+5)(n)*x(n)
        ple(P5, y_ptr, 5, EVP1);                      // y(m+5) = mem[y(m+5)]
        padd(P30, P30, P22, EC0);                     // acc += A(m+6)(n)*x(n)
        ple(P6, y_ptr, 6, EVP1);                      // y(m+6) = mem[y(m+6)]
        padd(P31, P31, P23, EC0);                     // acc += A(m+7)(n)*x(n)
        ple(P7, y_ptr, 7, EVP1);                      // y(m+7) = mem[y(m+7)]
        a_ptr = a_prev_ptr + 8*a_row_offset;
        pld(P9, beta_ptr, 0, EFP0);                   // beta
        x_ptr = (uintptr_t)x;
        pld(P8, alpha_ptr, 0, EFP0);                  // alpha
        pmul(P0, P0, P9, EC0);                        // y(m+0) = y(m+0)*beta
        pmul(P1, P1, P9, EC0);                        // y(m+1) = y(m+1)*beta
        pmul(P2, P2, P9, EC0);                        // y(m+2) = y(m+2)*beta
        pmul(P3, P3, P9, EC0);                        // y(m+3) = y(m+3)*beta
        pmul(P24, P24, P8, EC0);                      // acc(0) = acc(0)*alpha
        pmul(P25, P25, P8, EC0);                      // acc(1) = acc(1)*alpha
        pmul(P26, P26, P8, EC0);                      // acc(2) = acc(2)*alpha
        pmul(P27, P27, P8, EC0);                      // acc(3) = acc(3)*alpha
        padd(P0, P0, P24, EC0);                       // y(m+0) = beta*y(m+0) + acc(0)
        pmul(P28, P28, P8, EC0);                      // acc(4) = acc(4)*alpha
        padd(P1, P1, P25, EC0);                       // y(m+1) = beta*y(m+1) + acc(1)
        pmul(P29, P29, P8, EC0);                      // acc(5) = acc(5)*alpha
        padd(P2, P2, P26, EC0);                       // y(m+2) = beta*y(m+2) + acc(2)
        pmul(P30, P30, P8, EC0);                      // acc(6) = acc(6)*alpha
        padd(P3, P3, P27, EC0);                       // y(m+3) = beta*y(m+3) + acc(3)
        pmul(P31, P31, P8, EC0);                      // acc(7) = acc(7)*alpha
        pse(P0, y_ptr, 0, EVP1);                      // mem[y(m+0)] = y(m+0)
        pse(P1, y_ptr, 1, EVP1);                      // mem[y(m+1)] = y(m+1)
        pse(P2, y_ptr, 2, EVP1);                      // mem[y(m+2)] = y(m+2)
        pse(P3, y_ptr, 3, EVP1);                      // mem[y(m+3)] = y(m+3)
        pmul(P4, P4, P9, EC0);                        // y(m+4) = y(m+4)*beta
        pmul(P5, P5, P9, EC0);                        // y(m+5) = y(m+5)*beta
        pmul(P6, P6, P9, EC0);                        // y(m+6) = y(m+6)*beta
        pmul(P7, P7, P9, EC0);                        // y(m+7) = y(m+7)*beta
        padd(P4, P4, P28, EC0);                       // y(m+4) = beta*y(m+4) + acc(4)
        padd(P5, P5, P29, EC0);                       // y(m+5) = beta*y(m+5) + acc(5)
        padd(P6, P6, P30, EC0);                       // y(m+6) = beta*y(m+6) + acc(6)
        padd(P7, P7, P31, EC0);                       // y(m+7) = beta*y(m+7) + acc(7)
        pse(P4, y_ptr, 4, EVP1);                      // mem[y(m+4)] = y(m+4)
        pse(P5, y_ptr, 5, EVP1);                      // mem[y(m+5)] = y(m+5)
        pse(P6, y_ptr, 6, EVP1);                      // mem[y(m+6)] = y(m+6)
        pse(P7, y_ptr, 7, EVP1);                      // mem[y(m+7)] = y(m+7)
        y_ptr += 8*y_bytes;
    }

    //  process the remaining rows (not multiple of 8)
    if ((m % 8) != 0) {
        dvgemv_rows_chunk_1(m - i*8, n,
                            alpha,
                            (void*)a_ptr, lda,
                            (void*)x_ptr, x_bytes,
                            beta,
                            (void*)y_ptr, y_bytes);
    }
}


static void __dvgemv_cols_chunk_8_trans(int m,
                                        const double *a, int lda,
                                        const void *x, int x_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t x_ptr = (uintptr_t)x;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_offset = a_bytes*lda;

    ple(P8, x_ptr, 0, EVP0);            // x(m)
    ple(P9, x_ptr, 1, EVP0);            // x(m+1)
    x_ptr += 2*x_bytes;
    pld(P0, a_ptr, 0, EFP0);            // A(n)(m+0)
    pld(P1, a_ptr, 1, EFP0);            // A(n)(m+1)
    pmul(P16, P0, P8, EC0);             // A(n)(m+0)*x(m)
    pld(P2, a_ptr, 2, EFP0);            // A(n)(m+2)
    pld(P3, a_ptr, 3, EFP0);            // A(n)(m+3)
    pmul(P17, P1, P8, EC0);             // A(n)(m+1)*x(m)
    pld(P4, a_ptr, 4, EFP0);            // A(n)(m+4)
    pld(P5, a_ptr, 5, EFP0);            // A(n)(m+5)
    pmul(P18, P2, P8, EC0);             // A(n)(m+2)*x(m)
    pld(P6, a_ptr, 6, EFP0);            // A(n)(m+6)
    pld(P7, a_ptr, 7, EFP0);            // A(n)(m+7)
    pmul(P19, P3, P8, EC0);             // A(n)(m+3)*x(m)
    a_ptr += a_row_offset;
    for (int i = 0; i < (m - 1); i++) {
        pld(P0, a_ptr, 0, EFP0);        // A(n+1)(m+0)
        padd(P24, P24, P16, EC0);       // acc(0) += A(n)(m+0)*x(m)
        pmul(P20, P4, P8, EC0);         // A(n)(m+4)*x(m)
        pld(P1, a_ptr, 1, EFP0);        // A(n+1)(m+1)
        padd(P25, P25, P17, EC0);       // acc(1) += A(n)(m+1)*x(m)
        pmul(P21, P5, P8, EC0);         // A(n)(m+5)*x(m)
        pld(P2, a_ptr, 2, EFP0);        // A(n+1)(m+2)
        padd(P26, P26, P18, EC0);       // acc(2) += A(n)(m+2)*x(m)
        pmul(P22, P6, P8, EC0);         // A(n)(m+6)*x(m)
        pld(P3, a_ptr, 3, EFP0);        // A(n+1)(m+3)
        padd(P27, P27, P19, EC0);       // acc(3) += A(n)(m+3)*x(m)
        pmul(P23, P7, P8, EC0);         // A(n)(m+7)*x(m)
        pmv_p_p(P8, P9);                // x(m)
        ple(P9, x_ptr, 0, EVP0);        // x(m+1)
        x_ptr += x_bytes;
        pld(P4, a_ptr, 4, EFP0);        // A(n+1)(m+4)
        padd(P28, P28, P20, EC0);       // acc(4) += A(n+4)(m)*x(m)
        pmul(P16, P0, P8, EC0);         // A(n)(m+0)*x(m)
        pld(P5, a_ptr, 5, EFP0);        // A(n+1)(m+5)
        padd(P29, P29, P21, EC0);       // acc(5) += A(n+5)(m)*x(m)
        pmul(P17, P1, P8, EC0);         // A(n)(m+1)*x(m)
        pld(P6, a_ptr, 6, EFP0);        // A(n+1)(m+6)
        padd(P30, P30, P22, EC0);       // acc(6) += A(n)(m+6)*x(m)
        pmul(P18, P2, P8, EC0);         // A(n)(m+2)*x(m)
        pld(P7, a_ptr, 7, EFP0);        // A(n+1)(m+7)
        padd(P31, P31, P23, EC0);       // acc(7) += A(n)(m+7)*x(m)
        pmul(P19, P3, P8, EC0);         // A(n)(m+3)*x(m)
        a_ptr += a_row_offset;
    }
    padd(P24, P24, P16, EC0);           // acc(0) += A(n)(m+0)*x(m)
    pmul(P20, P4, P8, EC0);             // A(n+4)(m)*x(m)
    padd(P25, P25, P17, EC0);           // acc(1) += A(n)(m+1)*x(m)
    pmul(P21, P5, P8, EC0);             // A(n+5)(m)*x(m)
    padd(P26, P26, P18, EC0);           // acc(2) += A(n)(m+2)*x(m)
    pmul(P22, P6, P8, EC0);             // A(n+6)(m)*x(m)
    padd(P27, P27, P19, EC0);           // acc(3) += A(n)(m+3)*x(m)
    pmul(P23, P7, P8, EC0);             // A(n+7)(m)*x(m)
    padd(P28, P28, P20, EC0);           // acc(4) += A(n)(m+4)*x(m)
    padd(P29, P29, P21, EC0);           // acc(5) += A(n)(m+5)*x(m)
    padd(P30, P30, P22, EC0);           // acc(6) += A(n)(m+6)*x(m)
    padd(P31, P31, P23, EC0);           // acc(7) += A(n)(m+7)*x(m)
}


static void __dvgemv_cols_chunk_1_trans(int m,
                                        const double *a, int lda,
                                        const void *x, int x_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t x_ptr = (uintptr_t)x;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_offset = a_bytes*lda;

    ple(P8, x_ptr, 0, EVP0);            // x(m)
    x_ptr += x_bytes;
    pld(P0, a_ptr, 0, EFP0);            // A(n)(m)
    a_ptr += a_row_offset;
    pmul(P16, P0, P8, EC0);             // A(n)(m)*x(m)
    for (int i = 0; i < (m - 1); i++) {
        ple(P8, x_ptr, 0, EVP0);        // x(m+1)
        x_ptr += x_bytes;
        pld(P0, a_ptr, 0, EFP0);        // A(n+1)(m)
        a_ptr += a_row_offset;
        padd(P24, P24, P16, EC0);       // acc(0) += A(n)(m)*x(m)
        pmul(P16, P0, P8, EC0);         // A(n+1)(m)*x(m+1)
    }
    padd(P24, P24, P16, EC0);           // acc(0) += A(n)(m)*x(m)
}


static inline void __dvgemv_y_write_accumulators_alpha1_beta0(
        int n,
        void *y, int y_bytes)
{
    uintptr_t y_ptr = (uintptr_t)y;

    pse(P24, y_ptr, 0, EVP1);           // mem[y(n+0)] = y(n+0)
    if (n == 8) {
        pse(P25, y_ptr, 1, EVP1);           // mem[y(n+1)] = y(n+1)
        pse(P26, y_ptr, 2, EVP1);           // mem[y(n+2)] = y(n+2)
        pse(P27, y_ptr, 3, EVP1);           // mem[y(n+3)] = y(n+3)
        pse(P28, y_ptr, 4, EVP1);           // mem[y(n+4)] = y(n+4)
        pse(P29, y_ptr, 5, EVP1);           // mem[y(n+5)] = y(n+5)
        pse(P30, y_ptr, 6, EVP1);           // mem[y(n+6)] = y(n+6)
        pse(P31, y_ptr, 7, EVP1);           // mem[y(n+7)] = y(n+7)
    }
}


static inline void __dvgemv_y_write_accumulators(
        int n,
        const double alpha,
        const double beta,
        void       *y, int y_bytes)
{
    uintptr_t y_ptr = (uintptr_t)y;

    if ((alpha == 1.0) && (beta == 0.0)) {
        return __dvgemv_y_write_accumulators_alpha1_beta0(n, y, y_bytes);
    }

    pcvt_d_p(P8, dtoraw(alpha));            // alpha
    pcvt_d_p(P9, dtoraw(beta));             // beta

    if (n == 8) {
        pmul(P24, P24, P8, EC0);            // acc(0) = alpha*acc(0)
        ple(P0, y_ptr, 0, EVP1);            // y(n+0)
        pmul(P25, P25, P8, EC0);            // acc(1) = alpha*acc(1)
        ple(P1, y_ptr, 1, EVP1);            // y(n+1)
        pmul(P26, P26, P8, EC0);            // acc(2) = alpha*acc(2)
        ple(P2, y_ptr, 2, EVP1);            // y(n+2)
        pmul(P27, P27, P8, EC0);            // acc(3) = alpha*acc(3)
        ple(P3, y_ptr, 3, EVP1);            // y(n+3)
        pmul(P28, P28, P8, EC0);            // acc(4) = alpha*acc(4)
        ple(P4, y_ptr, 4, EVP1);            // y(n+4)
        pmul(P29, P29, P8, EC0);            // acc(5) = alpha*acc(5)
        ple(P5, y_ptr, 5, EVP1);            // y(n+5)
        pmul(P30, P30, P8, EC0);            // acc(6) = alpha*acc(6)
        ple(P6, y_ptr, 6, EVP1);            // y(n+6)
        pmul(P31, P31, P8, EC0);            // acc(7) = alpha*acc(7)
        ple(P7, y_ptr, 7, EVP1);            // y(n+7)
        pmul(P0, P0, P9, EC0);              // y(n+0) = beta*y(n+0)
        pmul(P1, P1, P9, EC0);              // y(n+1) = beta*y(n+1)
        pmul(P2, P2, P9, EC0);              // y(n+2) = beta*y(n+2)
        pmul(P3, P3, P9, EC0);              // y(n+3) = beta*y(n+3)
        padd(P0, P0, P24, EC0);             // y(n+0) = beta*y(n+0) + acc(0)
        padd(P1, P1, P25, EC0);             // y(n+1) = beta*y(n+1) + acc(1)
        pmul(P4, P4, P9, EC0);              // y(n+4) = beta*y(n+4)
        padd(P2, P2, P26, EC0);             // y(n+2) = beta*y(n+2) + acc(2)
        padd(P3, P3, P27, EC0);             // y(n+3) = beta*y(n+3) + acc(3)
        pse(P0, y_ptr, 0, EVP1);            // mem[y(n+0)] = y(n+0)
        pse(P1, y_ptr, 1, EVP1);            // mem[y(n+1)] = y(n+1)
        pmul(P5, P5, P9, EC0);              // y(n+5) = beta*y(n+5)
        padd(P4, P4, P28, EC0);             // y(n+4) = beta*y(n+4) + acc(4)
        padd(P5, P5, P29, EC0);             // y(n+5) = beta*y(n+5) + acc(5)
        pse(P2, y_ptr, 2, EVP1);            // mem[y(n+2)] = y(n+2)
        pse(P3, y_ptr, 3, EVP1);            // mem[y(n+3)] = y(n+3)
        pmul(P6, P6, P9, EC0);              // y(n+6) = beta*y(n+6)
        pse(P4, y_ptr, 4, EVP1);            // mem[y(n+4)] = y(n+4)
        pse(P5, y_ptr, 5, EVP1);            // mem[y(n+5)] = y(n+5)
        padd(P6, P6, P30, EC0);             // y(n+6) = beta*y(n+6) + acc(6)
        pmul(P7, P7, P9, EC0);              // y(n+7) = beta*y(n+7)
        pse(P6, y_ptr, 6, EVP1);            // mem[y(n+6)] = y(n+6)
        padd(P7, P7, P31, EC0);             // y(n+7) = beta*y(n+7) + acc
        pse(P7, y_ptr, 7, EVP1);            // mem[y(n+7)] = y(n+7)
    }

    if (n == 1) {
        pmul(P24, P24, P8, EC0);            // acc(0) = alpha*acc(0)
        ple(P0, y_ptr, 0, EVP1);            // y(n+0)
        pmul(P0, P0, P9, EC0);              // y(n+0) = beta*y(n+0)
        padd(P0, P0, P24, EC0);             // y(n+0) = beta*y(n+0) + acc(0)
        pse(P0, y_ptr, 0, EVP1);            // mem[y(n+0)] = y(n+0)
    }
}


inline int __bytes_to_cachelines(int bytes)
{
    return (bytes + BSP_CONFIG_DCACHE_LINE_BYTES - 1)/BSP_CONFIG_DCACHE_LINE_BYTES;
}

static void dvgemv_cols_chunk_8_trans_prefetch(
        int m, int n,
        const double alpha,
        const double *a, int lda,
        const void *x, int x_bytes,
        const double beta,
        void       *y, int y_bytes,
        char enable_prefetch)
{
    const int ROWS_PER_BLOCK = 32;
    const int COLS_PER_BLOCK = 8;
    const int COLS_NBLOCKS = m/ROWS_PER_BLOCK;
    const int ROWS_NBLOCKS = n/COLS_PER_BLOCK;
    const int REMAINING_ROWS = m - COLS_NBLOCKS*ROWS_PER_BLOCK;
    const int REMAINING_COLS = n - ROWS_NBLOCKS*COLS_PER_BLOCK;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_block_offset = a_bytes*lda*ROWS_PER_BLOCK;
    const uintptr_t a_col_block_offset = a_bytes*COLS_PER_BLOCK;
    const uintptr_t x_row_block_offset = x_bytes*ROWS_PER_BLOCK;
    const uintptr_t y_col_block_offset = y_bytes*COLS_PER_BLOCK;

    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;

#if VBLAS_ENABLE_HWPF
    const int HWPF_ENGINE_ID = 0;
    const int HWPF_ENGINE_REM_ID = 1;
    hwpf_engine_params_t hwpf_param;
    hwpf_param.stride = __bytes_to_cachelines(lda*sizeof(double));
    hwpf_param.nblocks = ROWS_PER_BLOCK;
    hwpf_param.nlines = __bytes_to_cachelines(COLS_PER_BLOCK*sizeof(double));

    hwpf_engine_params_t hwpf_param_rem;
    hwpf_param_rem.stride = __bytes_to_cachelines(lda*sizeof(double));
    hwpf_param_rem.nblocks = REMAINING_ROWS;
    hwpf_param_rem.nlines = __bytes_to_cachelines(COLS_PER_BLOCK*sizeof(double));

    hwpf_engine_throttle_t hwpf_throttle;
    hwpf_throttle.nwait = 8;
    hwpf_throttle.ninflight = 0;

#ifdef VBLAS_PREFETCH_X
    const int HWPF_ENGINE_X_ID = 2;
    const int HWPF_ENGINE_X_REM_ID = 3;
    hwpf_engine_params_t hwpf_param_x;
    hwpf_param_x.stride = 1;
    hwpf_param_x.nblocks = __bytes_to_cachelines(ROWS_PER_BLOCK*x_bytes);
    hwpf_param_x.nlines = 1;

    hwpf_engine_params_t hwpf_param_x_rem;
    hwpf_param_x_rem.stride = 1;
    hwpf_param_x_rem.nblocks = __bytes_to_cachelines(REMAINING_ROWS*x_bytes);
    hwpf_param_x_rem.nlines = 1;

    hwpf_engine_throttle_t hwpf_throttle_x;
    hwpf_throttle_x.nwait = COLS_PER_BLOCK;
    hwpf_throttle_x.ninflight = COLS_PER_BLOCK;
#endif

    int enable_hwpf =
            hwpf_param.stride
            && hwpf_param.nblocks
            && hwpf_param.nlines
            && enable_prefetch;

    if (enable_hwpf) {
#if 0
        static int trace_done = 0;
        if ( trace_done == 0 )  { printf("Transpose prefetch enabled.\n"); trace_done=1; 
        printf("Enabling HW prefetch: stride=%ld / nblocks=%ld / nlines=%ld\n",
                hwpf_param.stride,
                hwpf_param.nblocks,
                hwpf_param.nlines);
        printf("Matrix base address : %lx\n", a_ptr);
        printf("Matrix parameters: m=%d - n=%d - lda:%d\n", m, n, lda);

        }
#endif
        hwpf_set_params(HWPF_ENGINE_ID, &hwpf_param);
        hwpf_set_throttle(HWPF_ENGINE_ID, &hwpf_throttle);
#ifdef VBLAS_PREFETCH_X
        hwpf_set_params(HWPF_ENGINE_X_ID, &hwpf_param_x);
        hwpf_set_throttle(HWPF_ENGINE_X_ID, &hwpf_throttle_x);
#endif
        if (REMAINING_ROWS) {
            hwpf_set_params(HWPF_ENGINE_REM_ID, &hwpf_param_rem);
            hwpf_set_throttle(HWPF_ENGINE_REM_ID, &hwpf_throttle);
#ifdef VBLAS_PREFETCH_X
            hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
            hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle_x);
#endif
        }
    }
#if 0
    else {
        static int trace2_done = 0;
        if ( trace2_done == 0 )  { printf("Transpose prefetch disabled.\n"); trace2_done=1;  }
    }
#endif
#endif

    for (int j = 0; j < ROWS_NBLOCKS; j++) {
        pcvt_d_p(P24, 0);                       // acc(0) = 0
        pcvt_d_p(P25, 0);                       // acc(1) = 0
        pcvt_d_p(P26, 0);                       // acc(2) = 0
        pcvt_d_p(P27, 0);                       // acc(3) = 0
        pcvt_d_p(P28, 0);                       // acc(4) = 0
        pcvt_d_p(P29, 0);                       // acc(5) = 0
        pcvt_d_p(P30, 0);                       // acc(6) = 0
        pcvt_d_p(P31, 0);                       // acc(7) = 0

        //  Process all but the last block of the block column
        for (int i = 0; i < COLS_NBLOCKS-1; i++) {
#if VBLAS_ENABLE_HWPF
            //  Prefetch next block in the same block column
            if (enable_hwpf) {
                const uintptr_t _addr = a_ptr + a_row_block_offset;
                hwpf_wait(HWPF_ENGINE_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_ID, _addr, HWPF_ENABLE);

#ifdef VBLAS_PREFETCH_X
                const uintptr_t _addr_x = x_ptr + x_row_block_offset;
                hwpf_wait(HWPF_ENGINE_X_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
#endif
            }
#endif
            //  Process the current block
            __dvgemv_cols_chunk_8_trans(ROWS_PER_BLOCK,
                                        (double*)a_ptr, lda,
                                        (void*)x_ptr, x_bytes);

            a_ptr += a_row_block_offset;
            x_ptr += x_row_block_offset;
        }

#if VBLAS_ENABLE_HWPF
        if (enable_hwpf) {
            if (!REMAINING_ROWS) {
                //  Prefetch next block in the next block column
                const uintptr_t _addr = (uintptr_t)a + a_col_block_offset*(j + 1);
                hwpf_wait(HWPF_ENGINE_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_ID, _addr, HWPF_ENABLE);

#ifdef VBLAS_PREFETCH_X
                const uintptr_t _addr_x = (uintptr_t)x;
                hwpf_wait(HWPF_ENGINE_X_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
#endif
            } else {
                const uintptr_t _addr = a_ptr + a_row_block_offset;
                hwpf_wait(HWPF_ENGINE_REM_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_REM_ID, _addr, HWPF_ENABLE);

#ifdef VBLAS_PREFETCH_X
                const uintptr_t _addr_x = x_ptr + x_row_block_offset;
                hwpf_wait(HWPF_ENGINE_X_REM_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
#endif
            }
        }
#endif
        //  Process the last block of the block column
        __dvgemv_cols_chunk_8_trans(ROWS_PER_BLOCK,
                                    (double*)a_ptr, lda,
                                    (void*)x_ptr, x_bytes);

        a_ptr += a_row_block_offset;
        x_ptr += x_row_block_offset;

        //  Process the remaining rows in the block column
        if (REMAINING_ROWS) {
#if VBLAS_ENABLE_HWPF
            // Prefetch next block in the next block column
            if (enable_hwpf) {
                const uintptr_t _addr = (uintptr_t)a + a_col_block_offset*(j + 1);
                hwpf_wait(HWPF_ENGINE_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_ID, _addr, HWPF_ENABLE);

#ifdef VBLAS_PREFETCH_X
                const uintptr_t _addr_x = (uintptr_t)x;
                hwpf_wait(HWPF_ENGINE_X_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
#endif
            }
#endif

            __dvgemv_cols_chunk_8_trans(REMAINING_ROWS,
                                        (double*)a_ptr, lda,
                                        (void*)x_ptr, x_bytes);
        }

        //  Write into y
        __dvgemv_y_write_accumulators(COLS_PER_BLOCK, alpha, beta,
                (void*)y_ptr, y_bytes);

        a_ptr  = (intptr_t)a + a_col_block_offset*(j + 1);
        x_ptr = (uintptr_t)x;
        y_ptr += y_col_block_offset;
    }

    //  Process the remaining cols
    for (int j = 0; j < REMAINING_COLS; j++) {
        pcvt_d_p(P24, 0);                       // acc(0) = 0

        __dvgemv_cols_chunk_1_trans(m,
                                    (double*)a_ptr, lda,
                                    (void*)x_ptr, x_bytes);

        a_ptr += a_bytes;

        //  Write into y
        __dvgemv_y_write_accumulators(1, alpha, beta, (void*)y_ptr, y_bytes);
        y_ptr += y_bytes;
    }
}

void dvgemv_cols_chunk_8_trans(int m, int n,
                                      const double alpha,
                                      const double *a, int lda,
                                      const void *x, int x_bytes,
                                      const double beta,
                                      void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_offset = a_bytes*lda;

    int j;
    for (j = 0; j < n/8; j++) {
        a_prev_ptr = a_ptr;
        ple(P8, x_ptr, 0, EVP0);            // x(m)
        x_ptr += x_bytes;
        pld(P0, a_ptr, 0, EFP0);            // A(n)(m+0)
        pcvt_d_p(P24, 0);                   // acc(0) = 0
        pld(P1, a_ptr, 1, EFP0);            // A(n)(m+1)
        pcvt_d_p(P25, 0);                   // acc(1) = 0
        pmul(P16, P0, P8, EC0);             // A(n)(m+0)*x(m)
        pld(P2, a_ptr, 2, EFP0);            // A(n)(m+2)
        pcvt_d_p(P26, 0);                   // acc(2) = 0
        pld(P3, a_ptr, 3, EFP0);            // A(n)(m+3)
        pcvt_d_p(P27, 0);                   // acc(3) = 0
        pmul(P17, P1, P8, EC0);             // A(n)(m+1)*x(m)
        pld(P4, a_ptr, 4, EFP0);            // A(n)(m+4)
        pcvt_d_p(P28, 0);                   // acc(4) = 0
        pld(P5, a_ptr, 5, EFP0);            // A(n)(m+5)
        pcvt_d_p(P29, 0);                   // acc(5) = 0
        pmul(P18, P2, P8, EC0);             // A(n)(m+2)*x(m)
        pld(P6, a_ptr, 6, EFP0);            // A(n)(m+6)
        pcvt_d_p(P30, 0);                   // acc(6) = 0
        pld(P7, a_ptr, 7, EFP0);            // A(n)(m+7)
        pcvt_d_p(P31, 0);                   // acc(7) = 0
        pmul(P19, P3, P8, EC0);             // A(n)(m+3)*x(m)

        a_ptr += a_row_offset;
        for (int i = 0; i < (m - 1); i++) {
            pld(P0, a_ptr, 0, EFP0);        // A(n+1)(m+0)
            padd(P24, P24, P16, EC0);       // acc(0) += A(n)(m+0)*x(m)
            pmul(P20, P4, P8, EC0);         // A(n)(m+4)*x(m)
            pld(P1, a_ptr, 1, EFP0);        // A(n+1)(m+1)
            padd(P25, P25, P17, EC0);       // acc(1) += A(n)(m+1)*x(m)
            pmul(P21, P5, P8, EC0);         // A(n)(m+5)*x(m)
            pld(P2, a_ptr, 2, EFP0);        // A(n+1)(m+2)
            padd(P26, P26, P18, EC0);       // acc(2) += A(n)(m+2)*x(m)
            pmul(P22, P6, P8, EC0);         // A(n)(m+6)*x(m)
            pld(P3, a_ptr, 3, EFP0);        // A(n+1)(m+3)
            padd(P27, P27, P19, EC0);       // acc(3) += A(n)(m+3)*x(m)
            pmul(P23, P7, P8, EC0);         // A(n)(m+7)*x(m)
            pld(P4, a_ptr, 4, EFP0);        // A(n+1)(m+4)
            padd(P28, P28, P20, EC0);       // acc(4) += A(n+4)(m)*x(m)
            pmul(P16, P0, P8, EC0);         // A(n)(m+0)*x(m)
            pld(P5, a_ptr, 5, EFP0);        // A(n+1)(m+5)
            padd(P29, P29, P21, EC0);       // acc(5) += A(n+5)(m)*x(m)
            pmul(P17, P1, P8, EC0);         // A(n)(m+1)*x(m)
            pld(P6, a_ptr, 6, EFP0);        // A(n+1)(m+6)
            padd(P30, P30, P22, EC0);       // acc(6) += A(n)(m+6)*x(m)
            pmul(P18, P2, P8, EC0);         // A(n)(m+2)*x(m)
            pld(P7, a_ptr, 7, EFP0);        // A(n+1)(m+7)
            padd(P31, P31, P23, EC0);       // acc(7) += A(n)(m+7)*x(m)
            pmul(P19, P3, P8, EC0);         // A(n)(m+3)*x(m)
            a_ptr += a_row_offset;
            ple(P8, x_ptr, 0, EVP0);        // x(m)
            x_ptr += x_bytes;
        }
        padd(P24, P24, P16, EC0);           // acc(0) += A(n)(m+0)*x(m)
        pmul(P20, P4, P8, EC0);             // A(n+4)(m)*x(m)
        padd(P25, P25, P17, EC0);           // acc(1) += A(n)(m+1)*x(m)
        pmul(P21, P5, P8, EC0);             // A(n+5)(m)*x(m)
        padd(P26, P26, P18, EC0);           // acc(2) += A(n)(m+2)*x(m)
        pmul(P22, P6, P8, EC0);             // A(n+6)(m)*x(m)
        padd(P27, P27, P19, EC0);           // acc(3) += A(n)(m+3)*x(m)
        pmul(P23, P7, P8, EC0);             // A(n+7)(m)*x(m)
        padd(P28, P28, P20, EC0);           // acc(4) += A(n)(m+4)*x(m)
        pld(P8, alpha_ptr, 0, EFP0);        // alpha
        padd(P29, P29, P21, EC0);           // acc(5) += A(n)(m+5)*x(m)
        pld(P9, beta_ptr, 0, EFP0);         // beta
        padd(P30, P30, P22, EC0);           // acc(6) += A(n)(m+6)*x(m)
        a_ptr = a_prev_ptr + 8*a_bytes;
        padd(P31, P31, P23, EC0);           // acc(7) += A(n)(m+7)*x(m)
        x_ptr = (uintptr_t)x;
        pmul(P24, P24, P8, EC0);            // acc(0) = alpha*acc(0)
        ple(P0, y_ptr, 0, EVP1);            // y(n+0)
        pmul(P25, P25, P8, EC0);            // acc(1) = alpha*acc(1)
        ple(P1, y_ptr, 1, EVP1);            // y(n+1)
        pmul(P26, P26, P8, EC0);            // acc(2) = alpha*acc(2)
        ple(P2, y_ptr, 2, EVP1);            // y(n+2)
        pmul(P27, P27, P8, EC0);            // acc(3) = alpha*acc(3)
        ple(P3, y_ptr, 3, EVP1);            // y(n+3)
        pmul(P28, P28, P8, EC0);            // acc(4) = alpha*acc(4)
        ple(P4, y_ptr, 4, EVP1);            // y(n+4)
        pmul(P29, P29, P8, EC0);            // acc(5) = alpha*acc(5)
        ple(P5, y_ptr, 5, EVP1);            // y(n+5)
        pmul(P30, P30, P8, EC0);            // acc(6) = alpha*acc(6)
        ple(P6, y_ptr, 6, EVP1);            // y(n+6)
        pmul(P31, P31, P8, EC0);            // acc(7) = alpha*acc(7)
        ple(P7, y_ptr, 7, EVP1);            // y(n+7)
        pmul(P0, P0, P9, EC0);              // y(n+0) = beta*y(n+0)
        pmul(P1, P1, P9, EC0);              // y(n+1) = beta*y(n+1)
        pmul(P2, P2, P9, EC0);              // y(n+2) = beta*y(n+2)
        pmul(P3, P3, P9, EC0);              // y(n+3) = beta*y(n+3)
        padd(P0, P0, P24, EC0);             // y(n+0) = beta*y(n+0) + acc(0)
        padd(P1, P1, P25, EC0);             // y(n+1) = beta*y(n+1) + acc(1)
        pmul(P4, P4, P9, EC0);              // y(n+4) = beta*y(n+4)
        padd(P2, P2, P26, EC0);             // y(n+2) = beta*y(n+2) + acc(2)
        padd(P3, P3, P27, EC0);             // y(n+3) = beta*y(n+3) + acc(3)
        pse(P0, y_ptr, 0, EVP1);            // mem[y(n+0)] = y(n+0)
        pse(P1, y_ptr, 1, EVP1);            // mem[y(n+1)] = y(n+1)
        pmul(P5, P5, P9, EC0);              // y(n+5) = beta*y(n+5)
        padd(P4, P4, P28, EC0);             // y(n+4) = beta*y(n+4) + acc(4)
        padd(P5, P5, P29, EC0);             // y(n+5) = beta*y(n+5) + acc(5)
        pse(P2, y_ptr, 2, EVP1);            // mem[y(n+2)] = y(n+2)
        pse(P3, y_ptr, 3, EVP1);            // mem[y(n+3)] = y(n+3)
        pmul(P6, P6, P9, EC0);              // y(n+6) = beta*y(n+6)
        pse(P4, y_ptr, 4, EVP1);            // mem[y(n+4)] = y(n+4)
        pse(P5, y_ptr, 5, EVP1);            // mem[y(n+5)] = y(n+5)
        padd(P6, P6, P30, EC0);             // y(n+6) = beta*y(n+6) + acc(6)
        pmul(P7, P7, P9, EC0);              // y(n+7) = beta*y(n+7)
        pse(P6, y_ptr, 6, EVP1);            // mem[y(n+6)] = y(n+6)
        padd(P7, P7, P31, EC0);             // y(n+7) = beta*y(n+7) + acc
        pse(P7, y_ptr, 7, EVP1);            // mem[y(n+7)] = y(n+7)
        y_ptr += 8*y_bytes;
    }

    //  process the remaining rows (not multiple of 8)
    if ((n % 8) != 0) {
        dvgemv_cols_chunk_4_trans(m, n - j*8,
                                  alpha,
                                  (void*)a_ptr, lda,
                                  (void*)x_ptr, x_bytes,
                                  beta,
                                  (void*)y_ptr, y_bytes);
    }
}

void dvgemv_cols_chunk_4_trans(int m, int n,
                                      const double alpha,
                                      const double *a, int lda,
                                      const void *x, int x_bytes,
                                      const double beta,
                                      void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    const int a_bytes = sizeof(double);
    const uintptr_t a_row_offset = a_bytes*lda;


    int j;
    for (j = 0; j < n/4; j++) {
        a_prev_ptr = a_ptr;
        ple(P8, x_ptr, 0, EVP0);        // x(m)
        x_ptr += x_bytes;
        pld(P0, a_ptr, 0, EFP0);        // A(m)(n+0)
        pcvt_d_p(P24, 0);               // acc(0) = 0
        pld(P1, a_ptr, 1, EFP0);        // A(m)(n+1)
        pcvt_d_p(P25, 0);               // acc(1) = 0
        pmul(P16, P0, P8, EC0);         // A(n)(m+0)*x(m)
        pld(P2, a_ptr, 2, EFP0);        // A(m)(n+2)
        pcvt_d_p(P26, 0);               // acc(2) = 0
        pld(P3, a_ptr, 3, EFP0);        // A(m)(n+3)
        pcvt_d_p(P27, 0);               // acc(3) = 0
        pmul(P17, P1, P8, EC0);         // A(n)(m+1)*x(m)
        a_ptr += a_row_offset;
        for (int i = 0; i < (m - 1); i++) {
            pld(P0, a_ptr, 0, EFP0);    // A(n+1)(m+0)
            padd(P24, P24, P16, EC0);   // acc(0) += A(n)(m+0)*x(m)
            pmul(P18, P2, P8, EC0);     // A(n)(m+2)*x(m)
            pld(P1, a_ptr, 1, EFP0);    // A(n+1)(m+1)
            padd(P25, P25, P17, EC0);   // acc(1) += A(n)(m+1)*x(m)
            pmul(P19, P3, P8, EC0);     // A(n)(m+3)*x(m)
            pld(P2, a_ptr, 2, EFP0);    // A(n+1)(m+2)
            padd(P26, P26, P18, EC0);   // acc(2) += A(n)(m+2)*x(m)
            pmul(P16, P0, P8, EC0);     // A(n)(m+0)*x(m)
            pld(P3, a_ptr, 3, EFP0);    // A(n+1)(m+3)
            padd(P27, P27, P19, EC0);   // acc(3) += A(n)(m+3)*x(m)
            pmul(P17, P1, P8, EC0);     // A(n)(m+1)*x(m)
            ple(P8, x_ptr, 0, EVP0);    // x(m+1)
            a_ptr += a_row_offset;
            x_ptr += x_bytes;
        }

        pld(P9, alpha_ptr, 0, EFP0);    // alpha
        ple(P0, y_ptr, 0, EVP1);        // y(n+0)
        padd(P24, P24, P16, EC0);       // acc(0) += A(n+0)(m)*x(m)
        pmul(P18, P2, P8, EC0);         // A(n)(m+2)*x(m)
        ple(P1, y_ptr, 1, EVP1);        // y(n+1)
        padd(P25, P25, P17, EC0);       // acc(1) += A(n)(m+1)*x(m)
        pmul(P19, P3, P8, EC0);         // A(n)(m+3)*x(m)
        ple(P2, y_ptr, 2, EVP1);        // y(n+2)
        padd(P26, P26, P18, EC0);       // acc(2) += A(n)(m+2)*x(m)
        pmul(P24, P24, P9, EC0);        // acc(0) = alpha*acc(0)
        ple(P3, y_ptr, 3, EVP1);        // y(n+3)
        padd(P27, P27, P19, EC0);       // acc(3) += A(n)(m+3)*x(m)
        pmul(P25, P25, P9, EC0);        // acc(1) = alpha*acc(1)
        pld(P10, beta_ptr, 0, EFP0);    // beta
        a_ptr = a_prev_ptr + 4*a_bytes;
        x_ptr = (uintptr_t)x;
        pmul(P26, P26, P9, EC0);        // acc(2) = alpha*acc(2)
        pmul(P27, P27, P9, EC0);        // acc(3) = alpha*acc(3)
        pmul(P0, P0, P10, EC0);         // y(n+0) = beta*y(n+0)
        pmul(P1, P1, P10, EC0);         // y(n+1) = beta*y(n+1)
        pmul(P2, P2, P10, EC0);         // y(n+2) = beta*y(n+2)
        padd(P0, P0, P24, EC0);         // y(n+0) = beta*y(n+0) + alpha*acc(0)
        pmul(P3, P3, P10, EC0);         // y(n+3) = beta*y(n+3)
        padd(P1, P1, P25, EC0);         // y(n+1) = beta*y(n+1) + alpha*acc(1)
        padd(P2, P2, P26, EC0);         // y(n+2) = beta*y(n+2) + alpha*acc(2)
        padd(P3, P3, P27, EC0);         // y(m+3) = beta*y(n+3) + alpha*acc(3)
        pse(P0, y_ptr, 0, EVP1);        // mem[y(n+0)] = y(n+0)
        pse(P1, y_ptr, 1, EVP1);        // mem[y(n+1)] = y(n+1)
        pse(P2, y_ptr, 2, EVP1);        // mem[y(n+2)] = y(n+2)
        pse(P3, y_ptr, 3, EVP1);        // mem[y(n+3)] = y(n+3)
        y_ptr += 4*y_bytes;
    }

    //  process the remaining rows (not multiple of 8)
    if ((n % 4) != 0) {
        dvgemv_cols_chunk_1_trans(m, n - j*4,
                                  alpha,
                                  (double*)a_ptr, lda,
                                  (void*)x_ptr, x_bytes,
                                  beta,
                                  (void*)y_ptr, y_bytes);
    }
}

void dvgemv_cols_chunk_1_trans(int m, int n,
                                      const double alpha,
                                      const double *a, int lda,
                                      const void *x, int x_bytes,
                                      const double beta,
                                      void       *y, int y_bytes)
{
    uintptr_t a_ptr = (uintptr_t)a;
    uintptr_t a_prev_ptr;
    uintptr_t x_ptr = (uintptr_t)x;
    uintptr_t y_ptr = (uintptr_t)y;
    uintptr_t alpha_ptr = (uintptr_t)&alpha;
    uintptr_t beta_ptr = (uintptr_t)&beta;
    int a_bytes = sizeof(double);
    int a_row_offset = a_bytes*lda;

    pcvt_d_p(P24, 0);                   // acc = 0
    pld(P8, alpha_ptr, 0, EFP0);        // alpha
    pld(P9, beta_ptr, 0, EFP0);         // beta
    for (int j = 0; j < n; j++) {
        a_prev_ptr = a_ptr;
        ple(P1, x_ptr, 0, EVP0);            // x(m)
        pld(P0, a_ptr, 0, EFP0);            // A(n)(m)
        a_ptr += a_row_offset;
        x_ptr += x_bytes;
        for (int i = 0; i < (m - 1); i++) {
            pmul(P16, P0, P1, EC0);         // A(n)(m)*x(m)
            pld(P0, a_ptr, 0, EFP0);        // A(n)(m+1)
            a_ptr += a_row_offset;
            ple(P1, x_ptr, 0, EVP0);        // x(m+1)
            x_ptr += x_bytes;
            padd(P24, P24, P16, EC0);       // acc += A(n)(m)*x(m)
        }
        pmul(P16, P0, P1, EC0);             // A(n)(m)*x(m)
        x_ptr = (uintptr_t)x;
        a_ptr = a_prev_ptr + a_bytes;
        padd(P24, P24, P16, EC0);           // acc += A(n)(m)*x(m)
        pmul(P24, P24, P8, EC0);            // acc = alpha*acc
        ple(P0, y_ptr, 0, EVP1);            // y
        pmul(P0, P0, P9, EC0);              // y = beta*y
        padd(P0, P0, P24, EC0);             // y = y + acc
        pcvt_d_p(P24, 0);                   // acc = 0
        pse(P0, y_ptr, 0, EVP1);            // mem(y) = y(n)
        y_ptr += y_bytes;
    }
}

void dvgemv(int precision, char trans, int m, int n,
            const double alpha,
            const double *a, int lda,
            const void *x, vpfloat_evp_t x_evp,
            const double beta,
            void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* save environment */
    const uint64_t old_ec0 = pger_ec(EC0);
    const uint64_t old_efp0 = pger_efp(EFP0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    /* set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);
    pser_evp(pack_evp(
                x_evp.bis, vpfloat_get_rm_mem(),
                x_evp.es, x_evp.stride), EVP0);
    pser_evp(pack_evp(
                y_evp.bis, vpfloat_get_rm_mem(),
                y_evp.es, y_evp.stride), EVP1);

    if (trans == 'N') {
        if(lda <= (1<<16)){
            pser_efp(pack_efp(lda, vpfloat_get_rm_mem()), EFP0);
            dvgemv_rows_stride_chunk_8(
                m, n,
                alpha,
                a, lda,
                x, VPFLOAT_SIZEOF(x_evp),
                beta,
                y, VPFLOAT_SIZEOF(y_evp));
        } else {
            dvgemv_rows_chunk_8(
                m, n,
                alpha,
                a, lda,
                x, VPFLOAT_SIZEOF(x_evp),
                beta,
                y, VPFLOAT_SIZEOF(y_evp));            
        }
    } else {
        //dvgemv_cols_chunk_8_trans(
        dvgemv_cols_chunk_8_trans_prefetch(
                m, n,
                alpha,
                a, lda,
                x, VPFLOAT_SIZEOF(x_evp),
                beta,
                y, VPFLOAT_SIZEOF(y_evp),
                enable_prefetch);
    }

    /* restore environment */
    pser_ec(old_ec0, EC0);
    pser_efp(old_efp0, EFP0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    VBLASPERFMONITOR_FUNCTION_END;
}
