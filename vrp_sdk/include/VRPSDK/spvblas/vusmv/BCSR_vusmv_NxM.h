/**
* Copyright 2023 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        BCSR_vusmv_NxM.h
 *  @author      Valentin Isaac--Chassande
 */

#ifndef _BCSR_VUSMV_NXM_H_
#define _BCSR_VUSMV_NXM_H_

#include "Matrix/BCSR.h"

#define FUNC_PTR(R, C) (BCSR_dvusmv_ ## R ## x ## C)

void BCSR_dvusmv_NxM_handler(int precision, char trans,
                             const double alpha,
                             const dmatBCSR_t a,
                             const void * x, int x_bytes,
                             const double beta,
                             void * y, int y_bytes,
                             char enable_prefetch);

typedef struct _args_t {
    int precision; char trans;
    double alpha;
    dmatBCSR_t a;
    const void * x; int x_bytes;
    double beta;
    void * y; int y_bytes;
} args_t;

typedef void (*BCSR_dvusmv_NxM) (args_t args);

void FUNC_PTR(1,1)(args_t args);
void FUNC_PTR(2,1)(args_t args);
void FUNC_PTR(3,1)(args_t args);
void FUNC_PTR(4,1)(args_t args);
void FUNC_PTR(5,1)(args_t args);
void FUNC_PTR(6,1)(args_t args);
void FUNC_PTR(7,1)(args_t args);
void FUNC_PTR(8,1)(args_t args);

#endif /* _BCSR_VUSMV_NXM_H_ */
