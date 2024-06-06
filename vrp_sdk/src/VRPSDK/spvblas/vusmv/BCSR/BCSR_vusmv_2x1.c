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
 *  @file        BCSR_dvusmv_2x1.c
 *  @author      Valentin Isaac--Chassande
 */

#include <stdint.h>
#include "VRPSDK/asm/vpfloat.h"
#include "VRPSDK/vutils.h"
#include "VRPSDK/spvblas/vusmv/BCSR_vusmv_NxM.h"

static inline void _BCSR_dvusmv_2x1(args_t args)
{/* !!!! c * (2*1) */
    int c = args.a->col_block_size;
    
    int nbr = args.a->num_block_rows;

    int * bptr = args.a->bptr;
    int * bind = args.a->bind;
    double * bval = args.a->bval;

    uintptr_t x_ptr;
    uintptr_t y_ptr = (uintptr_t) args.y;
    uintptr_t bval_ptr = (uintptr_t) bval;

    int brow_start;
    int brow_end;

    pcvt_d_p(P0, dtoraw(args.alpha));
    pcvt_d_p(P1, dtoraw(args.beta));

    for (int brow = 0; brow < nbr; brow++) {
        brow_start = bptr[brow];
        brow_end = bptr[brow+1];

        pcvt_d_p(P30, 0);                                        // acc0 = 0
        pcvt_d_p(P31, 0);                                        // acc1 = 0

        for (int i = brow_start; i < brow_end; i++) {
            x_ptr = (uintptr_t) args.x + bind[i] * args.x_bytes;

            ple(P29, x_ptr, 0, EVP0);                            // x0

            pld(P2 , bval_ptr+0*c*sizeof(double), 0, EFP0);      // first bval0
            pld(P3 , bval_ptr+1*c*sizeof(double), 0, EFP0);      // first bval1

            bval_ptr += sizeof(double);
            x_ptr += args.x_bytes;
            for (int j = 0; j < c-1; j++) {

                pmul(P4 , P2 , P29, EC0);                            // bval0 * x0
                pld(P2 , bval_ptr+0*c*sizeof(double), 0, EFP0);      // next bval0
                pmul(P5 , P3 , P29, EC0);                            // bval1 * x0
                pld(P3 , bval_ptr+1*c*sizeof(double), 0, EFP0);      // next bval1

                ple(P29, x_ptr, 0, EVP0);                            // next x0

                bval_ptr += sizeof(double);
                x_ptr += args.x_bytes;

                padd(P30, P30, P4 , EC0);                            // acc0 += bval0*x0
                padd(P31, P31, P5 , EC0);                            // acc1 += bval1*x0
            }

            pmul(P4 , P2 , P29, EC0);                            // bval0 * x0
            pmul(P5 , P3 , P29, EC0);                            // bval1 * x0

            padd(P30, P30, P4 , EC0);                            // acc0 += bval0*x0
            padd(P31, P31, P5 , EC0);                            // acc1 += bval1*x0

            bval_ptr += (2*c-c)*sizeof(double);
        }

        pmul(P30, P30, P0, EC0);                                 // acc0 *= alpha
        ple(P2 , y_ptr, 0, EVP1);                                // y0
        pmul(P31, P31, P0, EC0);                                 // acc1 *= alpha
        ple(P3 , y_ptr, 1, EVP1);                                // y1

        pmul(P4 , P2 , P1, EC0);                                 // y0 *= beta
        pmul(P5 , P3 , P1, EC0);                                 // y1 *= beta

        padd(P4 , P4 , P30, EC0);                                // y0 += acc0
        pcvt_d_p(P30, 0);                                        // acc0 = 0
        padd(P5 , P5 , P31, EC0);                                // y1 += acc1
        pcvt_d_p(P31, 0);                                        // acc1 = 0

        pse(P4 , y_ptr, 0, EVP1);                                // mem(y0) = y0
        pse(P5 , y_ptr, 1, EVP1);                                // mem(y1) = y1

        y_ptr += 2*args.y_bytes;
    }
}

static inline void _BCSR_dvusmv_trans_2x1(args_t args)
{/* !!!! c * (2*1) */
    int c = args.a->col_block_size;

    int nbr = args.a->num_block_rows;

    int * bptr = args.a->bptr;
    int * bind = args.a->bind;
    double * bval = args.a->bval;

    uintptr_t x_ptr = (uintptr_t) args.x;
    uintptr_t y_ptr = (uintptr_t) args.y;
    uintptr_t bval_ptr = (uintptr_t) bval;

    int brow_start;
    int brow_end;

    pcvt_d_p(P0, dtoraw(args.alpha));

    for (int brow = 0; brow < nbr; brow++) {
        brow_start = bptr[brow];
        brow_end = bptr[brow+1];

        ple(P24, x_ptr, 0, EVP0);                                    // x0
        ple(P25, x_ptr, 1, EVP0);                                    // x1

        pmul(P24, P24, P0 , EC0);                                    // x0 *= alpha
        pmul(P25, P25, P0 , EC0);                                    // x1 *= alpha

        x_ptr += 2*args.x_bytes;
        for (int i = brow_start; i < brow_end; i++) {
            y_ptr = (uintptr_t) args.y + bind[i] * args.y_bytes;

            ple(P23, y_ptr, 0, EVP1);                            // y0

            pld(P1, bval_ptr+0*c*sizeof(double), 0, EFP0);       // bval0
            pld(P2, bval_ptr+1*c*sizeof(double), 0, EFP0);       // bval1

            bval_ptr += sizeof(double);
            for (int j = 0; j < c-1; j++) {

                pmul(P9 , P1, P24, EC0);                             // bval0 * x0
                pld(P1, bval_ptr+0*c*sizeof(double), 0, EFP0);       // bval0
                pmul(P10, P2, P25, EC0);                             // bval1 * x1
                pld(P2, bval_ptr+1*c*sizeof(double), 0, EFP0);       // bval1

                padd(P23, P23, P9 , EC0);                            // y0 += bval0*x0
                padd(P23, P23, P10, EC0);                            // y0 += bval1*x1

                pse(P23, y_ptr, 0, EVP1);                            // mem(y0) = y0

                bval_ptr += sizeof(double);
                y_ptr += args.y_bytes;

                ple(P23, y_ptr, 0, EVP1);                            // y0

            }

            pmul(P9 , P1, P24, EC0);                             // bval0 * x0
            pmul(P10, P2, P25, EC0);                             // bval1 * x1

            padd(P23, P23, P9 , EC0);                            // y0 += bval0*x0
            padd(P23, P23, P10, EC0);                            // y0 += bval1*x1

            pse(P23, y_ptr, 0, EVP1);                            // mem(y0) = y0

            bval_ptr += (2*c-c)*sizeof(double);
        }
    }
}

void BCSR_dvusmv_2x1(args_t args)
{
    if (args.trans == 'N')
        _BCSR_dvusmv_2x1(args);
    else
        _BCSR_dvusmv_trans_2x1(args);
}
