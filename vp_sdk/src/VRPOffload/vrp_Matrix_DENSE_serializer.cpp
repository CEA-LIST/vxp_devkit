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

#include "VRPOffload/vrp_Matrix_DENSE_serializer.hpp"

#define ALIGN_ADDRESS_64Bits(__offset) __offset = ( ( ( __offset / 8 ) + 1  ) * 8)

using namespace VPFloatPackage::Offloading;

size_t VRP_Matrix_DENSE_serializer::getSize(matrix_t a_matrix, size_t a_offset) {
    dmatDENSE_t l_dense = (dmatDENSE_t)a_matrix->matrix->repr;
    size_t l_size = 0;
    
    if ( l_dense != NULL ) {
        // Add size for lead_dim field
        l_size += sizeof(uint64_t);

        // Add size for val_size field
        l_size += sizeof(uint64_t);

        // padding
        l_size = (((( a_offset + l_size ) + VRP_Matrix_DENSE_serializer::getAlignment() - 1) / VRP_Matrix_DENSE_serializer::getAlignment()) * VRP_Matrix_DENSE_serializer::getAlignment() ) - a_offset;

        // compute valu field size for real data tyme
        size_t l_val_size = sizeof(double) * a_matrix->m * a_matrix->lda;

        // If val type is COMPLEX increase the size accordingly
        if ( a_matrix->type_value == COMPLEX_VALUE ) {    
            l_val_size *= 2;
        }

        // // Add size for val field
        l_size += l_val_size;
    }

    return l_size;
}

size_t VRP_Matrix_DENSE_serializer::getAlignment() {
    return 64;
}

void VRP_Matrix_DENSE_serializer::flaten(matrix_t a_matrix, uint64_t * a_free_address, uint64_t a_buffer_start_address) {
    dmatDENSE_t l_dense = (dmatDENSE_t)a_matrix->matrix->repr;
    uint64_t l_free_address = *a_free_address;
    uint64_t l_val_size = 0;
    uint64_t l_alignment = VRP_Matrix_DENSE_serializer::getAlignment();
    uint64_t l_padding = 0;

    // Add size for val field
    if ( a_matrix->type_value == COMPLEX_VALUE ) {
        l_val_size = sizeof(double) * 2 * ( a_matrix->m * a_matrix->lda );
    } else {
        l_val_size = sizeof(double) * ( a_matrix->m * a_matrix->lda );
    }

    memcpy((void *)l_free_address, &(l_dense->lead_dim), sizeof(int));
    l_free_address += sizeof(uint64_t);

    // Add field for val_size
    memcpy((void *)l_free_address, &(l_val_size), sizeof(uint64_t));
    l_free_address += sizeof(uint64_t);

    // Align val adress on cache size
    l_padding = (( (l_free_address - a_buffer_start_address + l_alignment - 1) / l_alignment ) * l_alignment); 

    l_free_address = a_buffer_start_address + l_padding;

    // Processing val field.
    memcpy((void *)l_free_address, l_dense->val,  l_val_size);
    l_free_address +=  l_val_size;

    // Update the address provided by caller
    *a_free_address = l_free_address;
}

dmatDENSE_t VRP_Matrix_DENSE_serializer::fromBuffer(uint64_t * a_address, uint64_t a_buffer_start_address) {
    dmatDENSE_t l_dense  = (dmatDENSE_t) malloc(sizeof(_dmatDENSE_t));
    uint64_t l_next_address = *a_address;
    uint64_t l_val_size;
    uint64_t l_alignment = VRP_Matrix_DENSE_serializer::getAlignment();
    uint64_t l_padding = 0;

    if ( l_dense == NULL) {
        std::cout << "Fail allocating memory for oski_matDENSE_t structure." << std::endl;
        return NULL;
    }

    l_dense->lead_dim = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_val_size = *(uint64_t *)l_next_address;
    l_next_address += sizeof(uint64_t);

    // Align val adress on cache size
    l_padding = (( (l_next_address + l_alignment - 1) / l_alignment ) * l_alignment) - a_buffer_start_address; 

    l_next_address = a_buffer_start_address + l_padding; 

    // Processing val field.
    l_dense->val = (double *)l_next_address;
    l_next_address += l_val_size ;

    // Update the address provided by caller
    *a_address = l_next_address;

    return l_dense;
}

void VRP_Matrix_DENSE_serializer::print(matrix_t a_matrix) {
    dmatDENSE_t l_dense = (dmatDENSE_t)a_matrix->matrix->repr;
    if ( l_dense ==  NULL ) {
        printf("NULL\n");
        return;
    }

    std::cout << "lead_dim : "<< l_dense->lead_dim << std::endl;

    std::cout << "val : ";
    for (int l_val_index = 0 ; l_val_index < a_matrix->m * a_matrix->n; l_val_index++) {
        std::cout << l_dense->val[l_val_index] << " ";
    }
    std::cout << std::endl;

}

