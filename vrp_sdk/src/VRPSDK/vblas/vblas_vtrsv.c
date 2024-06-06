/**
* Copyright 2022 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        vblas_vtrsv.c
 *  @author      Valentin Isaac--Chassande
 */
#include "VRPSDK/vblas.h"
#include "VRPSDK/vblas_perfmonitor.h"

static void vtrsv_1(int prec, char uplo, char diag, int n,
                    const void * a, vpfloat_evp_t a_evp, int lda,
                    void * x, vpfloat_evp_t x_evp)
{
    int a_bytes = VPFLOAT_SIZEOF(a_evp);
    int x_bytes = VPFLOAT_SIZEOF(x_evp);

    uintptr_t a_jptr = (uintptr_t) a;
    uintptr_t a_iptr = (uintptr_t) a;
    uintptr_t x_jptr = (uintptr_t) x;
    uintptr_t x_iptr = (uintptr_t) x;

    if (uplo == 'U' || uplo == 'u') {               // Upper triangular matrix

        a_jptr = (uintptr_t) a + (n-1)*(lda+1)*a_bytes;
	x_jptr = (uintptr_t) x + (n-1)*x_bytes;
        for (int j = n-1; j >= 0; j--) {
            ple(P0, x_jptr, 0, EVP1);               // x(j)
            a_iptr = (uintptr_t) a + ((j-1)*lda+j)*a_bytes;
            x_iptr = (uintptr_t) x + (j-1)*x_bytes;
            for (int i = j-1; i >= 0; i--) {
                ple(P1, a_iptr, 0, EVP0);           // a(i,j)
                ple(P2, x_iptr, 0, EVP0);           // x(i)
                pmul(P1, P1, P2, EC0);              // x(i)*a(i,j)
                a_iptr -= lda*a_bytes;
                psub(P0, P0, P1, EC0);              // x(j) -= x(i)*a(i,j)
                x_iptr -= x_bytes;
            }
            pse(P0, x_jptr, 0, EVP1);               // mem(x(j)) = x(j)
            if (diag == 'N' || diag == 'n') {       // if not unit diag
                vdiv(prec,
                     (void *) x_jptr, x_evp,
                     (void *) a_jptr, a_evp,
                     (void *) x_jptr, x_evp);       // mem(x(j)) = x(j)/a(j,j)
                a_jptr -= (lda+1)*a_bytes;
            }
            x_jptr -= x_bytes;
	}

    } else {                                        // Lower triangular matrix

        pcvt_d_p(P0, 0);
        for (int j = 0; j < n; j++) {
            ple(P1, x_jptr, 0, EVP1);               // x(j)
            if (pcmp_neq(P1, P0)) {                 // if (x(j) != 0)
                if (diag == 'N' || diag == 'n') {   // if not unit diag
                    vdiv(prec,
                         (void *) x_jptr, x_evp,
                         (void *) a_jptr, a_evp,
                         (void *) x_jptr, x_evp);   // mem(x(j)) = x(j)/a(j,j)
                    a_jptr += (lda+1)*a_bytes;
                    ple(P1, x_jptr, 0, EVP1);       // x(j)
                }
                x_iptr = x_jptr + x_bytes;
		a_iptr = (uintptr_t) a + (j+lda*(j+1))*a_bytes;
                for (int i = j+1; i < n; i++) {
                    ple(P3, a_iptr, 0, EVP0);       // a(i,j)
                    ple(P2, x_iptr, 0, EVP1);       // x(i)
                    pmul(P3, P1, P3, EC0);          // x(j)*a(i,j)
                    a_iptr += lda*a_bytes;
                    psub(P2, P2, P3, EC0);          // x(i) -= x(j)*a(i,j)
                    pse(P2, x_iptr, 0, EVP1);       // mem(x(i)) = x(i)
                    x_iptr += x_bytes;
                }
            }
            x_jptr += x_bytes;
	}

    }
}

static void vtrsv_trans_1(int prec, char uplo, char diag, int n,
                          const void * a, vpfloat_evp_t a_evp, int lda,
                          void * x, vpfloat_evp_t x_evp)
{
    int a_bytes = VPFLOAT_SIZEOF(a_evp);
    int x_bytes = VPFLOAT_SIZEOF(x_evp);

    uintptr_t a_jptr = (uintptr_t) a;
    uintptr_t a_iptr = (uintptr_t) a;
    uintptr_t x_jptr = (uintptr_t) x;
    uintptr_t x_iptr = (uintptr_t) x;

    if (uplo == 'U' || uplo == 'u') {               // Upper triangular matrix

        for (int j = 0; j < n; j++) {
            ple(P0, x_jptr, 0, EVP1);               // x(j)
            a_iptr = (uintptr_t) a + j*a_bytes;
            for (int i = 0; i < j; i++) {
                ple(P1, a_iptr, 0, EVP0);           // a(i,j)
                ple(P2, x_iptr, 0, EVP0);           // x(i)
                pmul(P1, P1, P2, EC0);              // x(i)*a(i,j)
                a_iptr += lda*a_bytes;
                psub(P0, P0, P1, EC0);              // x(j) -= x(i)*a(i,j)
                x_iptr += x_bytes;
            }
            pse(P0, x_jptr, 0, EVP1);               // mem(x(j)) = x(j)
            if (diag == 'N' || diag == 'n') {       // if not unit diag
                vdiv(prec,
                     (void *) x_jptr, x_evp,
                     (void *) a_jptr, a_evp,
                     (void *) x_jptr, x_evp);       // mem(x(j)) = x(j)/a(j,j)
                a_jptr += (lda+1)*a_bytes;
            }
            x_iptr = (uintptr_t) x;
            x_jptr += x_bytes;
	}

    } else {                                        // Lower triangular matrix

        pcvt_d_p(P0, 0);
        a_jptr = (uintptr_t) a + (n-1)*(lda+1)*a_bytes;
	x_jptr = (uintptr_t) x + (n-1)*a_bytes;
        for (int j = n-1; j >= 0; j--) {
            ple(P1, x_jptr, 0, EVP1);               // x(j)
            if (pcmp_neq(P1, P0)) {                 // if (x(j) != 0)
                if (diag == 'N' || diag == 'n') {   // if not unit diag
                    vdiv(prec,
                         (void *) x_jptr, x_evp,
                         (void *) a_jptr, a_evp,
                         (void *) x_jptr, x_evp);   // mem(x(j)) = x(j)/a(j,j)
                    a_jptr -= (lda+1)*a_bytes;
                    ple(P1, x_jptr, 0, EVP1);       // x(j)
                }
                x_iptr = (uintptr_t) x + (n-1)*x_bytes;
		a_iptr = (uintptr_t) a + (j+lda*(n-1))*a_bytes;
                for (int i = n-1; i > j; i--) {
                    ple(P3, a_iptr, 0, EVP0);       // a(i,j)
                    ple(P2, x_iptr, 0, EVP1);       // x(i)
                    pmul(P3, P1, P3, EC0);          // x(j)*a(i,j)
                    a_iptr -= lda*a_bytes;
                    psub(P2, P2, P3, EC0);          // x(i) -= x(j)*a(i,j)
                    pse(P2, x_iptr, 0, EVP1);       // mem(x(i)) = x(i)
                    x_iptr -= x_bytes;
                }
            }
            x_jptr -= x_bytes;
	}

    }
}

void vtrsv(int precision, char uplo, char trans, char diag, int n,
           const void * a, vpfloat_evp_t a_evp, int lda,
           void * x, vpfloat_evp_t x_evp)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* Save environment */
    const uint64_t old_ec0  = pger_ec(EC0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    /* Set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);
    pser_evp(pack_evp(a_evp.bis, vpfloat_get_rm_mem(),
             a_evp.es, a_evp.stride), EVP0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(),
             x_evp.es, x_evp.stride), EVP1);

    if (trans == 'N' || trans == 'n') {
        vtrsv_1(precision, uplo, diag, n,
                a, a_evp, lda,
                x, x_evp);
    } else {
        vtrsv_trans_1(precision, uplo, diag, n,
                      a, a_evp, lda,
                      x, x_evp);
    }

    /* Restore environment*/
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    VBLASPERFMONITOR_FUNCTION_END;
}

void strsv(char uplo, char trans, char diag, int n,
           const float * a, int lda,
           float * x, int x_inc)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_FLOAT;
    x_env.stride = x_inc;

    vtrsv(32, uplo, trans, diag, n,
          (void *) a, VPFLOAT_EVP_FLOAT, lda,
          (void *) x, x_env);

        VBLASPERFMONITOR_FUNCTION_END;
}

void dtrsv(char uplo, char trans, char diag, int n,
           const double * a, int lda,
           double * x, int x_inc)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_DOUBLE;
    x_env.stride = x_inc;

    vtrsv(64, uplo, trans, diag, n,
          (void *) a, VPFLOAT_EVP_DOUBLE, lda,
          (void *) x, x_env);
    
        VBLASPERFMONITOR_FUNCTION_END;
}
