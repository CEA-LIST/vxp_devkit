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
 *  @file        CSR_vusmv_NxM.h
 *  @author      Jerome Fereyre
 */

#ifndef _CSR_VUSMV_NXM_H_
#define _CSR_VUSMV_NXM_H_

#include "Matrix/CSR.h"

void CSR_dvusmv_NxM_handler(int precision, char trans, int m, int n,
                             const double alpha,
                             const dmatCSR_t a,
                             const void * x, int x_bytes,
                             const double beta,
                             void * y, int y_bytes,
                             char enable_prefetch);

void CSR_dvusmv_NxM_alpha1_beta0(int precision, char trans, int m, int n,
                             const dmatCSR_t a,
                             const void * x, int x_bytes,
                             void * y, int y_bytes,
                             char enable_prefetch);

#endif /* _CSR_VUSMV_NXM_H_ */

