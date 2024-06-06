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
 *  @file        vblas_mt_vgemv.c
 *  @author      Cesar Fuguet-Tortolero
 */
#include <stdio.h>
#include "VRPSDK/vblas.h"
#include "VRPSDK/vblas_mt.h"
#include "common/fifobuf.h"
#include "VRPSDK/asm/vpfloat.h"
#include "VRPSDK/vblas_perfmonitor.h"

typedef struct vblas_mt_vgemv_args_s
{
    int precision;
    char trans;
    int m;
    int n;
    const void *alpha;
    vpfloat_evp_t alpha_evp;
    const void *a;
    vpfloat_evp_t a_evp;
    int lda;
    const void *x;
    vpfloat_evp_t x_evp;
    const void *beta;
    vpfloat_evp_t beta_evp;
    void *y;
    vpfloat_evp_t y_evp;
    char enable_prefetch;
} vblas_mt_vgemv_args_t;

void vgemv_mt_kernel(void *args);

static inline int __vgemv_mt_compute_njobs(vblas_mt_vgemv_config_t *config,
                                           int m, int n)
{
    int njobs = m / config->min_rows_per_job;
    return njobs > 0 ? njobs : 1;
}

static inline fifobuf_node_t* __vgemv_mt_create_job_node(
        int job_id,
        int precision, char trans, int m, int n,
        const void *alpha, vpfloat_evp_t alpha_evp,
        const void *a, vpfloat_evp_t a_evp, int lda,
        const void *x, vpfloat_evp_t x_evp,
        const void *beta, vpfloat_evp_t beta_evp,
        void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    vblas_mt_vgemv_args_t *args;
    args = (vblas_mt_vgemv_args_t*)malloc(sizeof(vblas_mt_vgemv_args_t));
    args->precision = precision;
    args->trans = trans;
    args->m = m;
    args->n = n;
    args->alpha = alpha;
    args->alpha_evp = alpha_evp;
    args->a = a;
    args->a_evp = a_evp;
    args->lda = lda;
    args->x = x;
    args->x_evp = x_evp;
    args->beta = beta;
    args->beta_evp = beta_evp;
    args->y = y;
    args->y_evp = y_evp;
    args->enable_prefetch = enable_prefetch;

    vblas_job_t *job = (vblas_job_t*)malloc(sizeof(vblas_job_t));
    job->id = job_id;
    job->routine = vgemv_mt_kernel;
    job->args = args;

    fifobuf_node_t *job_node;
    job_node = (fifobuf_node_t*)malloc(sizeof(fifobuf_node_t));
    job_node->data = (void*)job;

    return job_node;
}

void vgemv_mt_kernel(void *args)
{
    vblas_mt_vgemv_args_t *vgemv_args = (vblas_mt_vgemv_args_t*)args;

    //  invalidate the entire cache to ensure cache coherency
    //  this is very painful but necessary without hardware cache coherency
    cpu_dcache_invalidate();

    vgemv(vgemv_args->precision,
          vgemv_args->trans,
          vgemv_args->m,
          vgemv_args->n,
          vgemv_args->alpha,
          vgemv_args->alpha_evp,
          vgemv_args->a,
          vgemv_args->a_evp,
          vgemv_args->lda,
          vgemv_args->x,
          vgemv_args->x_evp,
          vgemv_args->beta,
          vgemv_args->beta_evp,
          vgemv_args->y,
          vgemv_args->y_evp,
          vgemv_args->enable_prefetch);

    free(vgemv_args);
    cpu_dfence();
}

void vgemv_mt(vblas_mt_vgemv_config_t *config,
              int precision, char trans, int m, int n,
              const void *alpha, vpfloat_evp_t alpha_evp,
              const void *a, vpfloat_evp_t a_evp, int lda,
              const void *x, vpfloat_evp_t x_evp,
              const void *beta, vpfloat_evp_t beta_evp,
              void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    //  create and push distributed jobs into the job queue
    fifobuf_node_t *job_node;
    int njobs = __vgemv_mt_compute_njobs(config, m, n);
    int rows_per_job = m / njobs;
    int remaining_rows = m;
    int i;

    VBLASPERFMONITOR_FUNCTION_BEGIN;
    
    for (i = 0; i < njobs; i++) {
        uintptr_t _a, _y;
        int _m, _n;

        if ( trans == 'N' ) {
            _a = (uintptr_t)a + (uintptr_t)i*rows_per_job*lda*VPFLOAT_SIZEOF(a_evp);
            _y = (uintptr_t)y + (uintptr_t)i*rows_per_job*VPFLOAT_SIZEOF(y_evp);
            _m = (i < (njobs - 1)) ? rows_per_job : remaining_rows;
            _n = n;
            remaining_rows -= _m;
        } else {
            _a = (uintptr_t)a + (uintptr_t)i*rows_per_job*VPFLOAT_SIZEOF(a_evp);
            _y = (uintptr_t)y + (uintptr_t)i*rows_per_job*VPFLOAT_SIZEOF(y_evp);
            _m = m;
            _n = (i < (njobs - 1)) ? rows_per_job : remaining_rows;
            remaining_rows -= _n;
        }
        job_node = __vgemv_mt_create_job_node(
                i, precision, trans,
                _m, _n,
                alpha, alpha_evp,
                (void*)_a, a_evp, lda,
                x, x_evp,
                beta, beta_evp,
                (void*)_y, y_evp, enable_prefetch);

        fifobuf_push(config->vblas_mt->jobs, job_node);

#if 0
        printf("DEBUG: push node (%p) in list (%p) ; "
               "node->next = %p ; node->prev = %p ; "
               "list->next = %p ; list->prev = %p\n",
                &job_node->member, &config->vblas_mt->jobs->list,
                job_node->member.next,
                job_node->member.prev,
                config->vblas_mt->jobs->list.next,
                config->vblas_mt->jobs->list.prev);
#endif

        atomic_fetch_add(&config->vblas_mt->pending_jobs, 1);
    }

    //  The director thread also takes jobs
    fifobuf_t *jobs = config->vblas_mt->jobs;
    atomic_int *pending_jobs = &config->vblas_mt->pending_jobs;
    while (!fifobuf_is_empty(jobs)) {
        //  pop a job and check that another thread did not already took it
        fifobuf_node_t *node = fifobuf_pop(jobs);
        if (node == NULL) continue;

        //  call the job routine with the given arguments
        vblas_job_t *job = (vblas_job_t*)node->data;
        cpu_dcache_invalidate_range((uintptr_t)job, sizeof(vblas_job_t));
        (job->routine)(job->args);

        //  signal that the job is terminated
        atomic_fetch_add(pending_jobs, -1);

        //  free job and node memory
        free(job);
        free(node);
    }

    //  the director thread wait for all jobs to be terminated
    int pending;
    do {
        pending = atomic_fetch_or(pending_jobs, 0);
        cpu_delay(1ULL >> 9);
    } while (pending > 0);

    //  invalidate the entire cache to ensure cache coherency
    //  this is very painful but necessary without hardware cache coherency
    cpu_dcache_invalidate();

    VBLASPERFMONITOR_FUNCTION_END;

}
