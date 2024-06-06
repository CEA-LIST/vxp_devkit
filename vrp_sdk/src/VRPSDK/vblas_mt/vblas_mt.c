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

#include "VRPSDK/vblas_mt.h"
#include "common/fifobuf.h"
#include "common/threads.h"

static int vblas_mt_thread_idle(void *args);

static inline vblas_job_t *__vblas_mt_get_job(fifobuf_node_t *node)
{
    vblas_job_t *job = (vblas_job_t*)node->data;
    cpu_dcache_invalidate_range((uintptr_t)job, sizeof(vblas_job_t));
    return job;
}

int vblas_mt_init(vblas_mt_config_t *config)
{
    int ngood;

    //  create job queue
    config->jobs = (fifobuf_t*)malloc(sizeof(fifobuf_t));
    fifobuf_init(config->jobs);
    atomic_store(&config->pending_jobs, 0);
    cpu_dfence();

    //  create threads
    config->threads = malloc(config->max_threads*sizeof(thread_t));
    ngood = 0;
    for (int i = 0; i < config->max_threads; i++) {
        int thread_status;
        config->threads[i].cpu_id = NULL;
        thread_status = thread_create(&config->threads[i],
                vblas_mt_thread_idle,
                (void*)config);

        //  stop if one of the threads cannot be created
        if (thread_status < 0) break;
        ngood++;
    }

    //  if not all requested threads were created, destroy the created ones and
    //  return an error
    if (ngood < config->max_threads) {
        for (int i = 0; i < ngood; i++) {
            thread_destroy(&config->threads[i]);
        }
    }
    // return the number of threads that were not created as requested (as a
    // negative number)
    return ngood - config->max_threads;
}

int vblas_mt_destroy(vblas_mt_config_t *config)
{
    int nbad;

    for (int i = 0; i < config->max_threads; i++) {
        vblas_job_t *job;
        fifobuf_node_t *fifo_node;
        job = (vblas_job_t*)malloc(sizeof(vblas_job_t));
        job->id = 0;
        job->routine = NULL;
        job->args = NULL;
        fifo_node = (fifobuf_node_t*)malloc(sizeof(fifobuf_node_t));
        fifo_node->data = (void*)job;
        fifobuf_push(config->jobs, fifo_node);
    }

    nbad = 0;
    for (int i = 0; i < config->max_threads; i++) {
        if (thread_join(&config->threads[i]) < 0) nbad++;
    }

    free(config->jobs);
    free(config->threads);
    return -nbad;
}

static int vblas_mt_thread_idle(void *args)
{
    vblas_mt_config_t *config = (vblas_mt_config_t*)args;

    cpu_dcache_invalidate_range((uintptr_t)config, sizeof(vblas_mt_config_t*));
    fifobuf_t *jobs = config->jobs;

    for (;;) {
        //  wait for a job to be available
        if(fifobuf_is_empty(jobs)) continue;

        //  pop a job and check that another thread did not already took it
        fifobuf_node_t *node = fifobuf_pop(jobs);
        if (node == NULL) continue;

#if 0
        printf("DEBUG: pop node (%p) in list (%p) ; "
               "node->next = %p ; node->prev = %p ; "
               "list->next = %p ; list->prev = %p\n",
                &node->member, &jobs->list,
                node->member.next, node->member.prev,
                jobs->list.next, jobs->list.prev);
#endif

        //  stop thread if job routine is NULL
        vblas_job_t *job = __vblas_mt_get_job(node);
        if (job->routine == NULL) {
            free(job);
            free(node);
            break;
        }

        //  call the job routine with the given arguments
        (job->routine)(job->args);

        //  free job and node memory
        free(job);
        free(node);

        //  signal that the job is terminated
        atomic_fetch_add(&config->pending_jobs, -1);
    }
    return THREAD_SUCCESS;
}
