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

#ifndef __VRP_OFFLOADING_HPP__
#define __VRP_OFFLOADING_HPP__

#include "VRPOffload/vrp_argument_array.hpp"

#define VRP_OFFLAD_ENVIRONMENT_VAR_NAME "VRP_OFFLOAD"

namespace VPFloatPackage::Offloading {

    int call_solver(const VRPArgumentArray & a_argument_array, const char * a_solver_bin_path);

};

#endif /* __VRP_OFFLOADING_HPP__ */