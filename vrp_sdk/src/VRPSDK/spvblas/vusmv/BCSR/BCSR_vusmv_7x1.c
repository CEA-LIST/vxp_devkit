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
 *  @file        BCSR_dvusmv_7x1.c
 *  @author      Valentin Isaac--Chassande
 */

#include <stdint.h>
#include "VRPSDK/asm/vpfloat.h"
#include "VRPSDK/vutils.h"
#include "VRPSDK/spvblas/vusmv/BCSR_vusmv_NxM.h"

static inline void _BCSR_dvusmv_7x1(args_t args)
{/* !!!! c * (7*1) */
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

        pcvt_d_p(P24, 0);                                        // acc0 = 0
        pcvt_d_p(P25, 0);                                        // acc1 = 0
        pcvt_d_p(P26, 0);                                        // acc2 = 0
        pcvt_d_p(P27, 0);                                        // acc3 = 0
        pcvt_d_p(P28, 0);                                        // acc4 = 0
        pcvt_d_p(P29, 0);                                        // acc5 = 0
        pcvt_d_p(P30, 0);                                        // acc6 = 0

        for (int i = brow_start; i < brow_end; i++) {
            x_ptr = (uintptr_t) args.x + bind[i] * args.x_bytes;

            ple(P23, x_ptr, 0, EVP0);                            // x0

            pld(P2 , bval_ptr+0*c*sizeof(double), 0, EFP0);      // first bval0
            pld(P3 , bval_ptr+1*c*sizeof(double), 0, EFP0);      // first bval1
            pld(P4 , bval_ptr+2*c*sizeof(double), 0, EFP0);      // first bval2
            pld(P5 , bval_ptr+3*c*sizeof(double), 0, EFP0);      // first bval3
            pld(P6 , bval_ptr+4*c*sizeof(double), 0, EFP0);      // first bval4
            pld(P7 , bval_ptr+5*c*sizeof(double), 0, EFP0);      // first bval5
            pld(P8 , bval_ptr+6*c*sizeof(double), 0, EFP0);      // first bval6

            bval_ptr += sizeof(double);
            x_ptr += args.x_bytes;
            for (int j = 0; j < c-1; j++) {

                pmul(P10, P2 , P23, EC0);                            // bval0 * x0
                pld(P2 , bval_ptr+0*c*sizeof(double), 0, EFP0);      // next bval0
                pmul(P11, P3 , P23, EC0);                            // bval1 * x0
                pld(P3 , bval_ptr+1*c*sizeof(double), 0, EFP0);      // next bval1
                pmul(P12, P4 , P23, EC0);                            // bval2 * x0
                pld(P4 , bval_ptr+2*c*sizeof(double), 0, EFP0);      // next bval2
                pmul(P13, P5 , P23, EC0);                            // bval3 * x0
                pld(P5 , bval_ptr+3*c*sizeof(double), 0, EFP0);      // next bval3
                pmul(P14, P6 , P23, EC0);                            // bval4 * x0
                pld(P6 , bval_ptr+4*c*sizeof(double), 0, EFP0);      // next bval4
                pmul(P15, P7 , P23, EC0);                            // bval5 * x0
                pld(P7 , bval_ptr+5*c*sizeof(double), 0, EFP0);      // next bval5
                pmul(P16, P8 , P23, EC0);                            // bval6 * x0
                pld(P8 , bval_ptr+6*c*sizeof(double), 0, EFP0);      // next bval6

                ple(P23, x_ptr, 0, EVP0);                            // next x0

                bval_ptr += sizeof(double);
                x_ptr += args.x_bytes;

                padd(P24, P24, P10, EC0);                            // acc0 += bval0*x0
                padd(P25, P25, P11, EC0);                            // acc1 += bval1*x0
                padd(P26, P26, P12, EC0);                            // acc2 += bval2*x0
                padd(P27, P27, P13, EC0);                            // acc3 += bval3*x0
                padd(P28, P28, P14, EC0);                            // acc4 += bval4*x0
                padd(P29, P29, P15, EC0);                            // acc5 += bval5*x0
                padd(P30, P30, P16, EC0);                            // acc6 += bval6*x0
            }

            pmul(P10, P2 , P23, EC0);                            // bval0 * x0
            pmul(P11, P3 , P23, EC0);                            // bval1 * x0
            pmul(P12, P4 , P23, EC0);                            // bval2 * x0
            pmul(P13, P5 , P23, EC0);                            // bval3 * x0
            pmul(P14, P6 , P23, EC0);                            // bval4 * x0
            pmul(P15, P7 , P23, EC0);                            // bval5 * x0
            pmul(P16, P8 , P23, EC0);                            // bval6 * x0

            padd(P24, P24, P10, EC0);                            // acc0 += bval0*x0
            padd(P25, P25, P11, EC0);                            // acc1 += bval1*x0
            padd(P26, P26, P12, EC0);                            // acc2 += bval2*x0
            padd(P27, P27, P13, EC0);                            // acc3 += bval3*x0
            padd(P28, P28, P14, EC0);                            // acc4 += bval4*x0
            padd(P29, P29, P15, EC0);                            // acc5 += bval5*x0
            padd(P30, P30, P16, EC0);                            // acc6 += bval6*x0

            bval_ptr += (7*c-c)*sizeof(double);
        }

        pmul(P24, P24, P0, EC0);                                 // acc0 *= alpha
        ple(P2 , y_ptr, 0, EVP1);                                // y0
        pmul(P25, P25, P0, EC0);                                 // acc1 *= alpha
        ple(P3 , y_ptr, 1, EVP1);                                // y1
        pmul(P26, P26, P0, EC0);                                 // acc2 *= alpha
        ple(P4 , y_ptr, 2, EVP1);                                // y2
        pmul(P27, P27, P0, EC0);                                 // acc3 *= alpha
        ple(P5 , y_ptr, 3, EVP1);                                // y3
        pmul(P28, P28, P0, EC0);                                 // acc4 *= alpha
        ple(P6 , y_ptr, 4, EVP1);                                // y4
        pmul(P29, P29, P0, EC0);                                 // acc5 *= alpha
        ple(P7 , y_ptr, 5, EVP1);                                // y5
        pmul(P30, P30, P0, EC0);                                 // acc5 *= alpha
        ple(P8 , y_ptr, 6, EVP1);                                // y6

        pmul(P10, P2 , P1, EC0);                                 // y0 *= beta
        pmul(P11, P3 , P1, EC0);                                 // y1 *= beta
        pmul(P12, P4 , P1, EC0);                                 // y2 *= beta
        pmul(P13, P5 , P1, EC0);                                 // y3 *= beta
        pmul(P14, P6 , P1, EC0);                                 // y4 *= beta
        pmul(P15, P7 , P1, EC0);                                 // y5 *= beta
        pmul(P16, P8 , P1, EC0);                                 // y6 *= beta

        padd(P10, P10, P24, EC0);                                // y0 += acc0
        pcvt_d_p(P24, 0);                                        // acc0 = 0
        padd(P11, P11, P25, EC0);                                // y1 += acc1
        pcvt_d_p(P25, 0);                                        // acc1 = 0
        padd(P12, P12, P26, EC0);                                // y2 += acc2
        pcvt_d_p(P26, 0);                                        // acc2 = 0
        padd(P13, P13, P27, EC0);                                // y3 += acc3
        pcvt_d_p(P27, 0);                                        // acc3 = 0
        padd(P14, P14, P28, EC0);                                // y4 += acc4
        pcvt_d_p(P28, 0);                                        // acc4 = 0
        padd(P15, P15, P29, EC0);                                // y5 += acc5
        pcvt_d_p(P29, 0);                                        // acc5 = 0
        padd(P16, P16, P30, EC0);                                // y6 += acc6
        pcvt_d_p(P30, 0);                                        // acc6 = 0

        pse(P10, y_ptr, 0, EVP1);                                // mem(y0) = y0
        pse(P11, y_ptr, 1, EVP1);                                // mem(y1) = y1
        pse(P12, y_ptr, 2, EVP1);                                // mem(y2) = y2
        pse(P13, y_ptr, 3, EVP1);                                // mem(y3) = y3
        pse(P14, y_ptr, 4, EVP1);                                // mem(y4) = y4
        pse(P15, y_ptr, 5, EVP1);                                // mem(y5) = y5
        pse(P16, y_ptr, 6, EVP1);                                // mem(y6) = y6

        y_ptr += 7*args.y_bytes;
    }
}

static inline void _BCSR_dvusmv_trans_7x1(args_t args)
{/* !!!! c * (7*1) */
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
        ple(P26, x_ptr, 2, EVP0);                                    // x2
        ple(P27, x_ptr, 3, EVP0);                                    // x3
        ple(P28, x_ptr, 4, EVP0);                                    // x4
        ple(P29, x_ptr, 5, EVP0);                                    // x5
        ple(P30, x_ptr, 6, EVP0);                                    // x6

        pmul(P24, P24, P0 , EC0);                                    // x0 *= alpha
        pmul(P25, P25, P0 , EC0);                                    // x1 *= alpha
        pmul(P26, P26, P0 , EC0);                                    // x2 *= alpha
        pmul(P27, P27, P0 , EC0);                                    // x3 *= alpha
        pmul(P28, P28, P0 , EC0);                                    // x4 *= alpha
        pmul(P29, P29, P0 , EC0);                                    // x5 *= alpha
        pmul(P30, P30, P0 , EC0);                                    // x6 *= alpha

        x_ptr += 7*args.x_bytes;
        for (int i = brow_start; i < brow_end; i++) {
            y_ptr = (uintptr_t) args.y + bind[i] * args.y_bytes;

            ple(P23, y_ptr, 0, EVP1);                            // y0

            pld(P1, bval_ptr+0*c*sizeof(double), 0, EFP0);       // bval0
            pld(P2, bval_ptr+1*c*sizeof(double), 0, EFP0);       // bval1
            pld(P3, bval_ptr+2*c*sizeof(double), 0, EFP0);       // bval2
            pld(P4, bval_ptr+3*c*sizeof(double), 0, EFP0);       // bval3
            pld(P5, bval_ptr+4*c*sizeof(double), 0, EFP0);       // bval4
            pld(P6, bval_ptr+5*c*sizeof(double), 0, EFP0);       // bval5
            pld(P7, bval_ptr+6*c*sizeof(double), 0, EFP0);       // bval6

            bval_ptr += sizeof(double);
            for (int j = 0; j < c-1; j++) {

                pmul(P9 , P1, P24, EC0);                             // bval0 * x0
                pld(P1, bval_ptr+0*c*sizeof(double), 0, EFP0);       // bval0
                pmul(P10, P2, P25, EC0);                             // bval1 * x1
                pld(P2, bval_ptr+1*c*sizeof(double), 0, EFP0);       // bval1
                pmul(P11, P3, P26, EC0);                             // bval2 * x2
                pld(P3, bval_ptr+2*c*sizeof(double), 0, EFP0);       // bval2
                pmul(P12, P4, P27, EC0);                             // bval3 * x3
                pld(P4, bval_ptr+3*c*sizeof(double), 0, EFP0);       // bval3
                pmul(P13, P5, P28, EC0);                             // bval4 * x4
                pld(P5, bval_ptr+4*c*sizeof(double), 0, EFP0);       // bval4
                pmul(P14, P6, P29, EC0);                             // bval5 * x5
                pld(P6, bval_ptr+5*c*sizeof(double), 0, EFP0);       // bval5
                pmul(P15, P7, P30, EC0);                             // bval6 * x6
                pld(P7, bval_ptr+6*c*sizeof(double), 0, EFP0);       // bval6

                padd(P23, P23, P9 , EC0);                            // y0 += bval0*x0
                padd(P23, P23, P10, EC0);                            // y0 += bval1*x1
                padd(P23, P23, P11, EC0);                            // y0 += bval2*x2
                padd(P23, P23, P12, EC0);                            // y0 += bval3*x3
                padd(P23, P23, P13, EC0);                            // y0 += bval4*x4
                padd(P23, P23, P14, EC0);                            // y0 += bval5*x5
                padd(P23, P23, P15, EC0);                            // y0 += bval6*x6

                pse(P23, y_ptr, 0, EVP1);                            // mem(y0) = y0

                bval_ptr += sizeof(double);
                y_ptr += args.y_bytes;

                ple(P23, y_ptr, 0, EVP1);                            // y0

            }

            pmul(P9 , P1, P24, EC0);                             // bval0 * x0
            pmul(P10, P2, P25, EC0);                             // bval1 * x1
            pmul(P11, P3, P26, EC0);                             // bval2 * x2
            pmul(P12, P4, P27, EC0);                             // bval3 * x3
            pmul(P13, P5, P28, EC0);                             // bval4 * x4
            pmul(P14, P6, P29, EC0);                             // bval5 * x5
            pmul(P15, P7, P30, EC0);                             // bval6 * x6

            padd(P23, P23, P9 , EC0);                            // y0 += bval0*x0
            padd(P23, P23, P10, EC0);                            // y0 += bval1*x1
            padd(P23, P23, P11, EC0);                            // y0 += bval2*x2
            padd(P23, P23, P12, EC0);                            // y0 += bval3*x3
            padd(P23, P23, P13, EC0);                            // y0 += bval4*x4
            padd(P23, P23, P14, EC0);                            // y0 += bval5*x5
            padd(P23, P23, P15, EC0);                            // y0 += bval6*x6

            pse(P23, y_ptr, 0, EVP1);                            // mem(y0) = y0

            bval_ptr += (7*c-c)*sizeof(double);
        }
    }
}

void BCSR_dvusmv_7x1(args_t args)
{
    if (args.trans == 'N')
        _BCSR_dvusmv_7x1(args);
    else
        _BCSR_dvusmv_trans_7x1(args);
}
