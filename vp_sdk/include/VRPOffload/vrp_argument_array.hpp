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
 * Authors       : Alexandre Hoffman, Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#ifndef __VRP_ARGUMENT_ARRAY_HPP__
#define __VRP_ARGUMENT_ARRAY_HPP__

#include <stdint.h>
#include <deque>
#include "VRPSDK/spvblas.h"
#include "VPSDK/VBLASConfig.hpp"
#include "vrp_argument.hpp"

using namespace VPFloatPackage::VBLAS;

namespace VPFloatPackage::Offloading {

    class VRPArgumentArray {
        private:
            std::deque<VRPArgument> m_arguments;

        public:
            void addArgument(double * a_double_address, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(uint64_t * a_uint64_t_address, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(uint16_t * a_uint16_t_address, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(int32_t * a_int32_t_address, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(void * a_buffer_address, uint64_t a_buffer_size, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(matrix_t a_matrix_descr, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(VBLASConfig* a_vblas_config, vrp_argument_direction_t a_direction=VRP_SOLVER_ARGUMENT_IN_OUT);
            void addArgument(VRPArgument a_argument);

            friend vrp_solver_argument_array_t * marshall(const VRPArgumentArray & a_argument_array);
            friend VRPArgumentArray unmarshall(uint64_t a_source_address);
    };
};

#endif /* __VRP_ARGUMENT_ARRAY_HPP__ */