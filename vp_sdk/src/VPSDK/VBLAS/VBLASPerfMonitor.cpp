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
#include <iostream>
#include <iomanip>

std::ostream& operator<<(std::ostream& a_out, function_stats_t * a_function_stat) {
    a_out << std::setw(20) << a_function_stat->m_cumulated_duration;
    a_out << " | " << std::setw(20) << a_function_stat->m_call_count;
    return a_out;
}

std::ostream& operator<<(std::ostream& a_out, perf_map_entry_t * a_map_entry) {
    a_out << std::setw(20) << a_map_entry->m_tag;
    a_out << " | " << std::setw(20) << a_map_entry->m_cpu_id;
    a_out << " | " << &(a_map_entry->m_data);
    return a_out;
}

std::ostream& operator<<(std::ostream& a_out, perf_map_t * a_perf_map) {
    a_out << " m_nb_used_entries: " << a_perf_map->m_nb_used_entries << std::endl;
    a_out << std::setw(20) << "function name" <<  " | " << std::setw(20) << "cpu_id" << " | " << std::setw(20) << "duration" << " | " << std::setw(20) << "#call" << std::endl;
    for ( int i = 0; i < a_perf_map->m_nb_used_entries; i++ ) {
        a_out << &(a_perf_map->m_map_entries[i]);
        a_out << std::endl;
    }
    if (strncmp(a_perf_map->m_map_entries[PERF_MAP_MAX_ELEMENTS].m_tag, "dummy", PERF_MAP_KEY_MAX_LENGTH) == 0) {
        a_out << &(a_perf_map->m_map_entries[PERF_MAP_MAX_ELEMENTS]);
        a_out << std::endl;
    }
    return a_out;
}
