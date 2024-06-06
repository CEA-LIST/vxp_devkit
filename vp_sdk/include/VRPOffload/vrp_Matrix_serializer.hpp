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

#ifndef __VRP_MATRIX_SERIALIZER_HPP__
#define __VRP_MATRIX_SERIALIZER_HPP__

#include "Matrix/matrix.h"
#include <stdlib.h>
#include <stdint.h>

namespace VPFloatPackage::Offloading::VRP_Matrix_serializer {
        void * serialize(matrix_t a_matrix);            
        matrix_t unserialize(uint64_t * a_address);
        void print(matrix_t a_matrix);

        size_t getSize(matrix_t a_matrix);
        size_t getAlignment(matrix_t a_matrix);
        matrix_t fromBuffer(uint64_t * a_address);
        void flaten(matrix_t a_matrix, uint64_t * a_free_address);
};

#endif /* __VRP_MATRIX_SERIALIZER_HPP__ */