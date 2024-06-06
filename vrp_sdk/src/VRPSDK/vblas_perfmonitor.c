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

#include "VRPSDK/vblas_perfmonitor.h"
#include "common/cpu.h"
#include <stdio.h>

int getEntryIndexLocked(const char * a_tag, int a_cpu_id);
function_stats_t * getFunctionStat(const char * a_tag, int a_cpu_id);

static perf_map_t g_perf_map;

void perfInit(void) {
    printf("VBLAS perf monitoring structure initialized!\n");
    ticket_mutex_init(&g_perf_map.m_lock);
    
    strncpy(g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_tag, "dummy", PERF_MAP_KEY_MAX_LENGTH);
    g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_cpu_id = 0;

    g_perf_map.m_nb_used_entries = 0;

    g_perf_map.m_initialized = 1;
}


perf_map_t * getPerfMap(void) {
    return &g_perf_map;
}

void function_startAt(const char * a_tag, int a_cpu_id, clock_t a_start_time) {
    function_stats_t * l_function_stat = getFunctionStat(a_tag, a_cpu_id);
    l_function_stat->m_start_time = a_start_time;
    l_function_stat->m_call_count++;
}

void function_stopAt(const char * a_tag, int a_cpu_id, clock_t a_stop_time) {
    function_stats_t * l_function_stat = getFunctionStat(a_tag, a_cpu_id);
    l_function_stat->m_cumulated_duration += a_stop_time - l_function_stat->m_start_time;
}

int addEntryLocked(const char * a_tag, int a_cpu_id) {
    int l_new_entry_index = g_perf_map.m_nb_used_entries;

    // Save the tag is the next available element
    strncpy(g_perf_map.m_map_entries[l_new_entry_index].m_tag, a_tag, PERF_MAP_KEY_MAX_LENGTH);
    g_perf_map.m_map_entries[l_new_entry_index].m_cpu_id = a_cpu_id;

    // printf("Add %s %d at index %d\n", g_perf_map.m_map_entries[l_new_entry_index].m_tag, a_cpu_id, l_new_entry_index);

    g_perf_map.m_nb_used_entries++;

    return l_new_entry_index;
}

function_stats_t * getFunctionStat(const char * a_tag, int a_cpu_id) {
    cpu_dcache_invalidate();

    ticket_mutex_lock(&(g_perf_map.m_lock));

    int l_entry_index = getEntryIndexLocked(a_tag, a_cpu_id);

    // Entry was found ?
    if (  l_entry_index > g_perf_map.m_nb_used_entries ) {
        // No !
        if ( g_perf_map.m_nb_used_entries < PERF_MAP_MAX_ELEMENTS ) {
            // There is some room left for a new entry. Add a new entry
            l_entry_index = addEntryLocked(a_tag, cpu_id());   
        }
    }

    ticket_mutex_unlock(&(g_perf_map.m_lock));

    return &(g_perf_map.m_map_entries[l_entry_index].m_data);
}

int getEntryIndexLocked(const char * a_tag, int a_cpu_id) {
    
    int l_searched_index = PERF_MAP_MAX_ELEMENTS;
    int i = 0;

    if ( g_perf_map.m_initialized == 0 ) {
        printf("VBLAS perf monitoring structure NOT initialized!\n");
    }
    
    int l_last_valid_index = g_perf_map.m_nb_used_entries - 1;

    if ( g_perf_map.m_nb_used_entries != 0 ) {

        while(  i <= l_last_valid_index &&
                ( strncmp(g_perf_map.m_map_entries[i].m_tag, a_tag, PERF_MAP_KEY_MAX_LENGTH) != 0 ||
                g_perf_map.m_map_entries[i].m_cpu_id != a_cpu_id ) ) {
            i++;
        }

        if ( i <= l_last_valid_index ) {
            l_searched_index = i;
        }
    }

    return l_searched_index;
}

void display_perfmonitor(void) {
    printf(" g_perf_map.m_nb_used_entries : %d\n", g_perf_map.m_nb_used_entries);
    for ( int i = 0; i < g_perf_map.m_nb_used_entries; i++ ) { 
        printf("%20.20s | %20d | %20ld | %20d \n",  g_perf_map.m_map_entries[i].m_tag, 
                                                    g_perf_map.m_map_entries[i].m_cpu_id,
                                                    g_perf_map.m_map_entries[i].m_data.m_cumulated_duration,
                                                    g_perf_map.m_map_entries[i].m_data.m_call_count);
    }
    if (strncmp(g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_tag, "dummy", PERF_MAP_KEY_MAX_LENGTH) == 0) {
        printf("%20.20s | %20d | %20ld | %20d \n",  g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_tag,
                                                    g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_cpu_id,
                                                    g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_data.m_cumulated_duration,
                                                    g_perf_map.m_map_entries[PERF_MAP_MAX_ELEMENTS].m_data.m_call_count);
    } 
}