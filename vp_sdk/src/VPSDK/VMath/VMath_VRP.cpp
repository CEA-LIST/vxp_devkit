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

#include "VPSDK/VMath.hpp"
#include "VRPSDK/vmath.h"

using namespace VPFloatPackage;

VPFloat VMath::vsqrt(const VPFloat & a) {
    void * l_result = ::vsqrt(	VPFloatComputingEnvironment::get_precision(),
                a.getData(), a.getEnvironment(),
                NULL, a.getEnvironment());

    return VPFloat(l_result, a.getEnvironment().es, a.getEnvironment().bis, a.getEnvironment().stride, true);
}