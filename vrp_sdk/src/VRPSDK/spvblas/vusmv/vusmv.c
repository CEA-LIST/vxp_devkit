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
 *  @file        vusmv.c
 *  @author      Valentin Isaac--Chassande
 */

#include "VRPSDK/spvblas.h"
#include "VRPSDK/vmath.h"
#include "VRPSDK/vutils.h"
#include "VRPSDK/spvblas/vusmv/vusmv.h"
#include "VRPSDK/vblas_perfmonitor.h"

void dvusmv(int precision, char trans,
            const double alpha,
            const matrix_t a,
            const void * x, vpfloat_evp_t x_evp,
            const double beta,
            void * y, vpfloat_evp_t y_evp,
            char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* Save environment */
    const uint64_t old_ec0 = pger_ec(EC0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_efp0 = pger_efp(EFP0);

    /* Set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(),
             x_evp.es, x_evp.stride), EVP0);
    pser_evp(pack_evp(y_evp.bis, vpfloat_get_rm_mem(),
             y_evp.es, y_evp.stride), EVP1);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);

    int x_bytes = VPFLOAT_SIZEOF(x_evp);
    int y_bytes = VPFLOAT_SIZEOF(y_evp);
    
    /* Compute the correspondant format dvusmv */
    switch (a->matrix->type_id) {
        case BCSR:
            BCSR_dvusmv_NxM_handler(precision, trans,
                                    alpha,
                                    (dmatBCSR_t) a->matrix->repr,
                                    x, x_bytes,
                                    beta,
                                    y, y_bytes, enable_prefetch);
            break;
        case CSR:
            CSR_dvusmv_NxM_handler( precision, trans, a->m, a->n,
                                    alpha,
                                    (dmatCSR_t) a->matrix->repr,
                                    x, x_bytes,
                                    beta,
                                    y, y_bytes, enable_prefetch);
            break;
    }

    /* Restore environment */
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_efp(old_efp0, EFP0);

    VBLASPERFMONITOR_FUNCTION_END;
}
