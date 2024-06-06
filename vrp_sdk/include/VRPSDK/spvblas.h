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
 *  @file        spvblas.h
 *  @author      Valentin Isaac--Chassande
 */

#ifndef _SPVBLAS_H_
#define _SPVBLAS_H_

#include <stdlib.h>
#include "VRPSDK/asm/vpfloat.h"
#include "Matrix/matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Sparse matrix-vector multiplication
 *
 *  y := alpha*a.x + beta*y
 *
 *  Matrix in double precision
 *  Vectors in variable precision
 *  Scalar in double precision
 */
void dvusmv(int precision, char trans,
            const double alpha,
            const matrix_t a,
            const void * x, vpfloat_evp_t x_evp,
            const double beta,
            void * y, vpfloat_evp_t y_evp,
            char enable_prefetch);

#ifdef __cplusplus
}
#endif

#endif /* _SPVBLAS_H_ */
