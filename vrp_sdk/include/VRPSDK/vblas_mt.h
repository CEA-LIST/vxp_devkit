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
 *  @file        vblas_mt.h
 *  @author      Cesar Fuguet-Tortolero
 */
#ifndef __VBLAS_MT_H__
#define __VBLAS_MT_H__

#ifdef __cplusplus
#include <atomic>
using namespace std;
extern "C" {
#else
#include <stdatomic.h>
#endif
#include "asm/vpfloat.h"

struct fifobuf_s;
struct thread_s;

/**
 *  Configuration of the VBLAS multi-thread server
 */
typedef struct vblas_mt_config_s
{
    int max_threads;
    struct fifobuf_s *jobs;
    atomic_int pending_jobs;
    struct thread_s *threads;
} vblas_mt_config_t;

/**
 *  Definition a VBLAS multi-thread job
 */
typedef struct vblas_job_s
{
    int  id;
    void (*routine)(void*);
    void *args;
} vblas_job_t;

/**
 *  Definition a VBLAS multi-thread VGEMV routine
 */
typedef struct vblas_mt_vgemv_config_s
{
    vblas_mt_config_t *vblas_mt;
    int min_rows_per_job;
} vblas_mt_vgemv_config_t;

/**
 *  @func   vblas_mt_init
 *  @brief  Initialize the server of threads for multi-threaded VBLAS routines
 *  @return Number of threads that were not created as requested (as a
 *          negative number)
 */
int vblas_mt_init(vblas_mt_config_t *config);

/**
 *  @func   vblas_mt_destroy
 *  @brief  Destroy the server of threads for multi-threaded VBLAS routines
 *  @return Number of threads that finished with an error (as a negative number)
 */
int vblas_mt_destroy(vblas_mt_config_t *config);

/**
 *  @func  vgemv_mt
 *  @brief Multi-threaded version of Matrix-Vector multiplication
 */
void vgemv_mt(vblas_mt_vgemv_config_t *config,
              int precision, char trans, int m, int n,
              const void *alpha, vpfloat_evp_t alpha_evp,
              const void *a, vpfloat_evp_t a_evp, int lda,
              const void *x, vpfloat_evp_t x_evp,
              const void *beta, vpfloat_evp_t beta_evp,
              void *y, vpfloat_evp_t y_evp, char enable_prefetch);

#ifdef __cplusplus
}
#endif

#endif /* __VBLAS_MT_H__ */
