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

#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__

#include <stdlib.h>
#include "Matrix/matrix.h"

namespace VPFloatPackage::Solver {

    int bicg(int precision, int transpose, int n, double * x, matrix_t A, matrix_t At, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);

    int precond_bicg(int precision, int transpose, int n, double * x, matrix_t A, matrix_t At, matrix_t iM, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);

    int bicgstab(int precision, int transpose, int n, double * x, matrix_t A, matrix_t At, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);

    int cg(int precision, int transpose, int n, double * x, matrix_t A, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);

    int precond_cg(int precision, int transpose, int n, double * x, matrix_t A, matrix_t iM, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);

    int qmr(int precision, int transpose, int n, double * x, matrix_t A, matrix_t At, double * b, double tolerance, uint16_t exponent_size = 7, int32_t stride_size = 1, char * log_buffer = NULL, uint64_t log_buffer_size = 0);
};

#endif /* __SOLVERS_HPP__ */