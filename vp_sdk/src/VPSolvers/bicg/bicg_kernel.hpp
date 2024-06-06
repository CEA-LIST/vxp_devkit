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

#ifndef __BICG_KERNEL_HPP__
#define __BICG_KERNEL_HPP__

#include "VPSDK/VPFloat.hpp"
#include "Matrix/matrix.h"

using namespace VPFloatPackage;

int bicg_vp(int precision, int transpose, int n, VPFloatArray & x,  matrix_t A, matrix_t At, VPFloatArray b, double tolerance, uint16_t exponent_size, int32_t stride_size);

#endif /*  __BICG_KERNEL_HPP__ */