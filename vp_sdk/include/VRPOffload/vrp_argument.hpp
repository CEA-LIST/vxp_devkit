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

#ifndef __VRP_ARGUMENT_HPP__
#define __VRP_ARGUMENT_HPP__

#include <stdint.h>
#include "vrp_ioctl.h"

namespace VPFloatPackage::Offloading {

    struct VRPArgument {
        VRPArgument(uint64_t a_size, uint64_t a_address, uint64_t a_alignment, vrp_argument_direction_t a_direction);
        VRPArgument(vrp_solver_argument_t  * a_vrp_solver_argument_address);
        VRPArgument();

        uint64_t m_size;
        uint64_t m_address;
        uint64_t m_alignment;
        vrp_argument_direction_t m_direction;
    };

};

#endif /* __VRP_ARGUMENT_HPP__ */