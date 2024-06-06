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
 * Authors       : Jerome Fereyre
 * Creation Date : October, 2023
 * Description   : 
 **/

#ifndef __VBLAS_PERFMONITOR_H__
#define __VBLAS_PERFMONITOR_H__

#include <time.h>
#include "common/ticket_mutex.h"
#include "common/cpu.h"

#define PERF_MAP_KEY_MAX_LENGTH 20
#define PERF_MAP_MAX_ELEMENTS   20

typedef struct function_stats {
    clock_t m_cumulated_duration;
    clock_t m_start_time;
    int     m_call_count;
} function_stats_t;

typedef struct perf_map_entry {
    char                m_tag[PERF_MAP_KEY_MAX_LENGTH];
    int                 m_cpu_id;
    function_stats_t    m_data;
} perf_map_entry_t;

typedef struct perf_map {
    int                 m_initialized;
    ticket_mutex_t      m_lock;
    atomic_int          m_nb_used_entries __cl_aligned__;
    perf_map_entry_t    m_map_entries[PERF_MAP_MAX_ELEMENTS+1];
} perf_map_t;

#ifdef __cplusplus
#include <iostream>
std::ostream& operator<<(std::ostream& a_out, perf_map_t * a_perf_map);

extern "C" {
#endif

void function_startAt(const char * a_tag, int a_cpu_id, clock_t a_start_time);
void function_stopAt(const char * a_tag, int a_cpu_id, clock_t a_stop_time);
void display_perfmonitor(void);
perf_map_t * getPerfMap(void);
void perfInit(void);


#ifdef __cplusplus
}
#endif

#if VBLAS_ENABLE_PERFMONITOR == 1

#define VBLASPERFMONITOR_BEGIN(function_)   function_startAt(function_, cpu_id(), clock());
#define VBLASPERFMONITOR_END(function_)     function_stopAt(function_, cpu_id(), clock()); 
#define VBLASPERFMONITOR_FUNCTION_BEGIN     VBLASPERFMONITOR_BEGIN(__FUNCTION__);
#define VBLASPERFMONITOR_FUNCTION_END       VBLASPERFMONITOR_END(__FUNCTION__); 
#define VBLASPERFMONITOR_INITIALIZE         perfInit();

#ifdef __cplusplus
#define VBLASPERFMONITOR_DISPLAY            std::cout << getPerfMap() << std::endl;
#else
#define VBLASPERFMONITOR_DISPLAY            display_perfmonitor();

#endif

#else 

#define VBLASPERFMONITOR_BEGIN
#define VBLASPERFMONITOR_END
#define VBLASPERFMONITOR_FUNCTION_BEGIN
#define VBLASPERFMONITOR_FUNCTION_END
#define VBLASPERFMONITOR_INITIALIZE
#define VBLASPERFMONITOR_DISPLAY            

#endif /* VBLAS_ENABLE_PERFMONITOR */



#endif /* __VBLAS_PERFMONITOR_H__ */