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

#include <stdlib.h>
#include <stdio.h>
#include "VRPOffload/vrp_argument.hpp"


VPFloatPackage::Offloading::VRPArgument::VRPArgument()  : m_size(0), m_alignment(0), m_address(0), m_direction(VRP_SOLVER_ARGUMENT_IN_OUT) {

}

VPFloatPackage::Offloading::VRPArgument::VRPArgument(uint64_t a_size, uint64_t a_address, uint64_t a_alignment, vrp_argument_direction_t a_direction) : m_size(a_size), m_alignment(a_alignment), m_address(a_address), m_direction(a_direction) {

}

VPFloatPackage::Offloading::VRPArgument::VRPArgument( vrp_solver_argument_t * a_vrp_solver_argument_address )
    : m_size(a_vrp_solver_argument_address->size)
    , m_alignment(a_vrp_solver_argument_address->alignment)
    , m_address((uint64_t)(a_vrp_solver_argument_address->address))
    , m_direction(a_vrp_solver_argument_address->direction) {
    
}