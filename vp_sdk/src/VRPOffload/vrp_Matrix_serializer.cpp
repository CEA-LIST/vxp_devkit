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
#include <string.h>
#include <iostream>

#include "VRPOffload/vrp_Matrix_serializer.hpp"
#include "VRPOffload/vrp_Matrix_CSR_serializer.hpp"
#include "VRPOffload/vrp_Matrix_BCSR_serializer.hpp"
#include "VRPOffload/vrp_Matrix_DENSE_serializer.hpp"

using namespace VPFloatPackage::Offloading;

size_t VRP_Matrix_serializer::getSize(matrix_t a_matrix) {
    size_t l_size = 0;

    // Sizeof of m, n, base_index, lda, format field
    l_size += sizeof(uint64_t) * 5;

    // Sizeof of type_matrix  field
    l_size += sizeof(uint64_t);

    // Sizeof of type_value  field
    l_size += sizeof(uint64_t);

    // Sizeof of matrix.type_id  field
    l_size += sizeof(uint64_t);

    switch(a_matrix->type_matrix) {
        case DENSE:
            l_size += VRP_Matrix_DENSE_serializer::getSize(a_matrix, l_size);
            break;
        case CSR:
            l_size += VRP_Matrix_CSR_serializer::getSize(a_matrix, l_size);
            break;
        case BCSR:
            l_size += VRP_Matrix_BCSR_serializer::getSize(a_matrix);
            break;
        default:
            std::cout << "Matrix type " << a_matrix->type_matrix << " not supported. Size is O" << std::endl;
            break; 
    }
    
    return l_size;
}

size_t VRP_Matrix_serializer::getAlignment(matrix_t a_matrix) {
        switch(a_matrix->type_matrix) {
        case DENSE:        
            return VRP_Matrix_DENSE_serializer::getAlignment();
            break;
        case CSR:        
            return VRP_Matrix_CSR_serializer::getAlignment();
            break;            
        default:
            return 0;
            break;
    }
}

void VRP_Matrix_serializer::flaten(matrix_t a_matrix, uint64_t * a_free_address) {
    uint64_t l_free_address = *a_free_address;
    uint64_t l_buffer_start_address = *a_free_address;

    memcpy((void *)l_free_address, &(a_matrix->m), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->n), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->base_index), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->lda), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->format), sizeof(char));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->type_matrix), sizeof(types_e));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->type_value), sizeof(types_value_e));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_matrix->matrix->type_id), sizeof(int));
    l_free_address += sizeof(uint64_t);

    switch(a_matrix->type_matrix) {
        case DENSE:        
            VRP_Matrix_DENSE_serializer::flaten(a_matrix, &l_free_address, l_buffer_start_address );
            break;
        case CSR:
            VRP_Matrix_CSR_serializer::flaten(a_matrix, &l_free_address, l_buffer_start_address );
            break;
        case BCSR:
            VRP_Matrix_BCSR_serializer::flaten(a_matrix, &l_free_address);
            break;
        default:
            std::cout << "Matrix type " << a_matrix->type_matrix << " not supported. flaten is partial." << std::endl;
            break; 
    }

    // Update the address provided by caller
    *a_free_address = l_free_address;
}

matrix_t VRP_Matrix_serializer::fromBuffer(uint64_t * a_address) {
    matrix_t l_matrix  = (_matrix_t *) malloc(sizeof(_matrix_t));
    uint64_t l_next_address = *a_address;
    uint64_t l_buffer_start_address = *a_address;

    if ( l_matrix == NULL ) {
        std::cout << "Fail allocating memory for _matrix_t structure." << std::endl;
        return NULL;
    }

    l_matrix->matrix  = (_oski_mat_t *) malloc(sizeof(_oski_mat_t));

    if ( l_matrix->matrix == NULL ) {
        std::cout << "Fail allocating memory for _oski_mat_t structure." << std::endl;
        free(l_matrix);
        return NULL;
    }

    l_matrix->m = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->n = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->base_index = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->lda = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->format = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->type_matrix = *(types_e *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->type_value = *(types_value_e *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_matrix->matrix->type_id = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    switch(l_matrix->type_matrix) {
        case DENSE:
            l_matrix->matrix->repr = (void *)VRP_Matrix_DENSE_serializer::fromBuffer(&l_next_address, l_buffer_start_address);
            break;
        case CSR:
            l_matrix->matrix->repr = (void *)VRP_Matrix_CSR_serializer::fromBuffer(&l_next_address, l_buffer_start_address);
            break;
        case BCSR:
            l_matrix->matrix->repr = (void *)VRP_Matrix_BCSR_serializer::fromBuffer(&l_next_address);
            break;
        default:
            std::cout << "Matrix type " << l_matrix->type_matrix << " not supported. fromBuffer is partial." << std::endl;
            break; 
    }

    // Update the address provided by caller
    *a_address = l_next_address;

    return l_matrix;
}

void * VRP_Matrix_serializer::serialize(matrix_t a_matrix) {
    size_t l_matrix_size = getSize(a_matrix);

    void * l_flat_matrix = malloc(l_matrix_size);

    uint64_t l_buffer_offset = (uint64_t)l_flat_matrix;

    flaten(a_matrix, &l_buffer_offset);

    return l_flat_matrix;
}

matrix_t VRP_Matrix_serializer::unserialize(uint64_t * a_address) {
    uint64_t l_offset_address = *a_address;
    
    matrix_t l_matrix = fromBuffer(&l_offset_address);

    *a_address = l_offset_address;

    return l_matrix;
}

void VRP_Matrix_serializer::print(matrix_t a_matrix) {
    if ( a_matrix ==  NULL ) {
        printf("NULL\n");
        return;
    }

    std::cout << "m : "<< a_matrix->m << std::endl;
    std::cout << "n : "<< a_matrix->m << std::endl;
    std::cout << "type_matrix : "<< a_matrix->type_matrix << std::endl;
    std::cout << "type_value : "<< a_matrix->type_value << std::endl;

    switch(a_matrix->type_matrix) {
        case DENSE:
            VRP_Matrix_DENSE_serializer::print(a_matrix);
            break;
        case CSR:
            VRP_Matrix_CSR_serializer::print(a_matrix);
            break;
        case BCSR:
            VRP_Matrix_BCSR_serializer::print(a_matrix);
            break;
        default:
            std::cout << "Matrix type : "<< a_matrix->type_matrix << " not supported." << std::endl;
    }
}

