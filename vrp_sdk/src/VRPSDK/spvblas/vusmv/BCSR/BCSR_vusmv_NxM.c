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
 *  @file        BCSR_vusmv_NxM.c
 *  @author      Valentin Isaac--Chassande
 */

#include <stdio.h>
#include <stddef.h>
#include "VRPSDK/spvblas/vusmv/BCSR_vusmv_NxM.h"
#include "VRPSDK/vblas_perfmonitor.h"

static BCSR_dvusmv_NxM BCSR_dvusmv_funcptr[8][8] =
{
    {
    &FUNC_PTR(1,1), &FUNC_PTR(1,1), &FUNC_PTR(1,1), &FUNC_PTR(1,1),
    &FUNC_PTR(1,1), &FUNC_PTR(1,1), &FUNC_PTR(1,1), &FUNC_PTR(1,1)
    },
    {
    &FUNC_PTR(2,1), &FUNC_PTR(2,1), &FUNC_PTR(2,1), &FUNC_PTR(2,1),
    &FUNC_PTR(2,1), &FUNC_PTR(2,1), &FUNC_PTR(2,1), &FUNC_PTR(2,1)
    },
    {
    &FUNC_PTR(3,1), &FUNC_PTR(3,1), &FUNC_PTR(3,1), &FUNC_PTR(3,1),
    &FUNC_PTR(3,1), &FUNC_PTR(3,1), &FUNC_PTR(3,1), &FUNC_PTR(3,1)
    },
    {
    &FUNC_PTR(4,1), &FUNC_PTR(4,1), &FUNC_PTR(4,1), &FUNC_PTR(4,1),
    &FUNC_PTR(4,1), &FUNC_PTR(4,1), &FUNC_PTR(4,1), &FUNC_PTR(4,1)
    },
    {
    &FUNC_PTR(5,1), &FUNC_PTR(5,1), &FUNC_PTR(5,1), &FUNC_PTR(5,1),
    &FUNC_PTR(5,1), &FUNC_PTR(5,1), &FUNC_PTR(5,1), &FUNC_PTR(5,1)
    },
    {
    &FUNC_PTR(6,1), &FUNC_PTR(6,1), &FUNC_PTR(6,1), &FUNC_PTR(6,1),
    &FUNC_PTR(6,1), &FUNC_PTR(6,1), &FUNC_PTR(6,1), &FUNC_PTR(6,1)
    },
    {
    &FUNC_PTR(7,1), &FUNC_PTR(7,1), &FUNC_PTR(7,1), &FUNC_PTR(7,1),
    &FUNC_PTR(7,1), &FUNC_PTR(7,1), &FUNC_PTR(7,1), &FUNC_PTR(7,1)
    },
    {
    &FUNC_PTR(8,1), &FUNC_PTR(8,1), &FUNC_PTR(8,1), &FUNC_PTR(8,1),
    &FUNC_PTR(8,1), &FUNC_PTR(8,1), &FUNC_PTR(8,1), &FUNC_PTR(8,1)
    }
};

void BCSR_dvusmv_NxM_handler(int precision, char trans,
                             const double alpha,
                             const dmatBCSR_t a,
                             const void * x, int x_bytes,
                             const double beta,
                             void * y, int y_bytes,
                             char enable_prefetch)
{
    int r = a->row_block_size;
    int c = a->col_block_size;

    BCSR_dvusmv_NxM funcptr = BCSR_dvusmv_funcptr[r-1][c-1];

    if (funcptr == 0) {
        printf("No implementation found\n");
        return;
    }

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    args_t args;
    args.precision = precision;
    args.trans     = trans;
    args.alpha     = alpha;
    args.a         = a;
    args.x         = x;
    args.x_bytes   = x_bytes;
    args.beta      = beta;
    args.y         = y;
    args.y_bytes   = y_bytes;

    (*funcptr) (args);

    if (a->num_rows_leftover && a->leftover != NULL) {

        if (trans == 'N')
            y += y_bytes * r * a->num_block_rows;
        else
            x += x_bytes * r * a->num_block_rows;

        BCSR_dvusmv_NxM_handler(precision, trans,
                                alpha,
                                (dmatBCSR_t) a->leftover,
                                x, x_bytes,
                                beta,
                                y, y_bytes,
                                enable_prefetch);

    }

    VBLASPERFMONITOR_FUNCTION_BEGIN;
}
