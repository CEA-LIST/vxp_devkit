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
 * Creation Date : August, 2023
 * Description   : 
 **/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include "VRPOffload/vrp_argument_array.hpp"
#include "VRPOffload/vrp_Matrix_serializer.hpp"
#include "vrp_ioctl.h"

using namespace VPFloatPackage::Offloading;

void VRPArgumentArray::addArgument(double * a_address, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(VRPArgument(sizeof(double), (uint64_t)a_address, 64, a_direction));
}

void VRPArgumentArray::addArgument(uint64_t * a_address, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(VRPArgument(sizeof(uint64_t), (uint64_t)a_address, 64, a_direction));
}

void VRPArgumentArray::addArgument(uint16_t * a_address, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(VRPArgument(sizeof(uint16_t), (uint64_t)a_address, 64, a_direction));
}

void VRPArgumentArray::addArgument(int32_t * a_address, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(VRPArgument(sizeof(int32_t), (uint64_t)a_address, 64, a_direction));
}

void VRPArgumentArray::addArgument(void * a_address, uint64_t a_buffer_size, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(VRPArgument(a_buffer_size, (uint64_t)a_address, 64, a_direction));
}

void VRPArgumentArray::addArgument(VRPArgument a_argument) {
    this->m_arguments.push_back(a_argument);
}

void VRPArgumentArray::addArgument(matrix_t a_matrix_descr, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(
        VRPArgument(VRP_Matrix_serializer::getSize(a_matrix_descr), 
                    (uint64_t)VRP_Matrix_serializer::serialize(a_matrix_descr), 
                    VRP_Matrix_serializer::getAlignment(a_matrix_descr),
                    a_direction
        )
    );
}

void VRPArgumentArray::addArgument(VBLASConfig* a_vblas_config, vrp_argument_direction_t a_direction) {
    this->m_arguments.push_back(
        VRPArgument(sizeof(VBLASConfig), 
                    (uint64_t)VBLAS_getConfig(), 
                    64,
                    a_direction
        )
    );
}

namespace VPFloatPackage::Offloading {

    vrp_solver_argument_array_t * marshall(const VRPArgumentArray & a_argument_array) {

        int l_argument_index = 0;
        vrp_solver_argument_array_t * l_argument_array = (vrp_solver_argument_array_t *)malloc(sizeof_vrp_solver_arguments(a_argument_array.m_arguments.size()));

        if ( l_argument_array == NULL ) {
            std::cout << "Fail allocating memory for argument array." << std::endl;
            return NULL;
        }

        l_argument_array->nb_arguments = a_argument_array.m_arguments.size();
        
        for (VRPArgument l_vrp_argument : a_argument_array.m_arguments ) {
            l_argument_array->arguments[l_argument_index].address = (void *)l_vrp_argument.m_address;
            l_argument_array->arguments[l_argument_index].size = l_vrp_argument.m_size;
            l_argument_array->arguments[l_argument_index].alignment = l_vrp_argument.m_alignment;
            l_argument_array->arguments[l_argument_index].direction = l_vrp_argument.m_direction;
            l_argument_index++;
        }

        return l_argument_array;
    }

    VRPArgumentArray unmarshall(uint64_t a_source_address) {
        VRPArgumentArray l_argument_array;
        uint64_t l_offset = sizeof(vrp_count_t);
        vrp_count_t l_argument_index = 0;

        vrp_count_t l_nb_arguments = *((vrp_count_t *)a_source_address);

        for ( l_argument_index = 0 ; l_argument_index < l_nb_arguments; l_argument_index++ ) {
            
            l_argument_array.addArgument( VRPArgument( (vrp_solver_argument_t *)(a_source_address + l_offset) ) );
            l_offset += sizeof(vrp_solver_argument_t);

        }

        return l_argument_array;
    }
}