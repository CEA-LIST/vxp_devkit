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
 * Authors       : Alexandre Hoffmann, Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#ifndef __QMR_KERNEL_HPP__
#define __QMR_KERNEL_HPP__

#include "VPSDK/VPFloat.hpp"
#include "Matrix/matrix.h"

using namespace VPFloatPackage;

int qmr_vp(               // Solves Ax = b where A is a non-symetric square matrix using the Quasi-minimal residual method from [1,2]    
	int precision,          // precision used by VPFloat (aka the size of the mantissa)
	int transpose,			// flag used to specify that matrices are transposed
	int n,                  // size of the matrix 
	VPFloatArray& x,        // an initial guess for the system Ax = b
	matrix_t A,             // A
	matrix_t At,            // Transpose of A
	VPFloatArray b,         // LHS
	double tolerance,       // stoping criterion for the method. The method will stop if |r_k|^2 < tolerance^2 
	uint16_t exponent_size, // size of the exponent of a VPFloat
	int32_t stride_size);   // mysterious variable which shall be 1
#endif /*  __QMR_KERNEL_HPP__ */

// [1] Freund, R. W., & Nachtigal, N. M. (1994). An implementation of the QMR method based on coupled two-term recurrences. SIAM Journal on Scientific Computing, 15(2), 313-337.
// [2] Barrett, R., Berry, M., Chan, T. F., Demmel, J., Donato, J., Dongarra, J., ... & Van der Vorst, H. (1994). Templates for the solution of linear systems: building blocks for iterative methods. Society for Industrial and Applied Mathematics.
