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

#include "VRPOffload/vrp_Matrix_CSR_serializer.hpp"

using namespace VPFloatPackage::Offloading;

#define ALIGN_ADDRESS_64Bits(__offset) __offset = ( ( ( __offset  + 7 ) / 8 ) * 8 )

/**
 * @brief 
 * 
 * @param l_csr 
 * @return size_t 
 */
size_t VRP_Matrix_CSR_serializer::getSize(matrix_t a_matrix, size_t a_offset) {
    dmatCSR_t l_csr = (dmatCSR_t)a_matrix->matrix->repr;
    size_t l_size = 0;
    size_t l_ptr_size = 0, l_ind_size = 0, l_val_size = 0;
    size_t l_alignment = VRP_Matrix_CSR_serializer::getAlignment();

    if ( l_csr != NULL ) {
        // Add size for base_index
        l_size += sizeof(uint64_t);

        // Add size for has_unit_diag_implicit
        l_size += sizeof(uint64_t);
        
        // Add size for has_sorted_indeces
        l_size += sizeof(uint64_t);
        
        // Add size for is_upper
        l_size += sizeof(uint64_t);
        
        // Add size for is_lower
        l_size += sizeof(uint64_t);
        
        // Processing ptr field.
        // Add size for ptr size
        l_size += sizeof(uint64_t);

        // padding
        l_size = ( ( l_size + l_alignment - 1 ) / l_alignment ) * l_alignment;

        // Add size for ptr field
        l_ptr_size = sizeof(int) * (a_matrix->m + 1);
        l_size += l_ptr_size;

        // padding
        l_size = ( ( l_size + 7 ) / 8 ) * 8;

        l_size += sizeof(uint64_t);

        // Add size for ind field
        l_ind_size = sizeof(int) * (l_csr->ptr[a_matrix->m] - l_csr->base_index);
        
        // padding
        l_size = ( ( l_size + l_alignment - 1 ) / l_alignment ) * l_alignment;
       
        // padding
        l_size += l_ind_size;

        // padding
        l_size = ( ( l_size + 7 ) / 8 ) * 8;

        // Add size for val field
        if ( a_matrix->type_value == COMPLEX_VALUE ) {
            l_val_size = sizeof(double) * 2 * (l_csr->ptr[a_matrix->m] - l_csr->base_index);
        } else{
            l_val_size = sizeof(double) * (l_csr->ptr[a_matrix->m] - l_csr->base_index);
        }

        l_size += sizeof(uint64_t);

        // padding
        l_size = ( ( l_size + l_alignment - 1 ) / l_alignment ) * l_alignment;

        // padding
        l_size += l_val_size;

        l_size = ( ( l_size + 7 ) / 8 ) * 8;

        // Add size for is_shared
        l_size += sizeof(uint64_t);
    }

    return l_size;
}

size_t VRP_Matrix_CSR_serializer::getAlignment() {
    return 64;
}

/**
 * @brief 
 * 
 * @param a_matrix 
 * @param a_address 
 */
void VRP_Matrix_CSR_serializer::flaten(matrix_t a_matrix, uint64_t * a_free_address, uint64_t a_buffer_start_address) {
    dmatCSR_t l_csr = (dmatCSR_t)a_matrix->matrix->repr;

    uint64_t l_free_address = *a_free_address;
    uint64_t l_ptr_size, l_ind_size, l_val_size;
    uint64_t l_alignment = VRP_Matrix_CSR_serializer::getAlignment();
    uint64_t l_padding = 0;

    memcpy((void *)l_free_address, &(l_csr->base_index), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(l_csr->has_unit_diag_implicit), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(l_csr->has_sorted_indices), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(l_csr->stored.is_upper), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(l_csr->stored.is_lower), sizeof(int));
    l_free_address += sizeof(uint64_t);

    // Processing ptr field.
    // 2 field are created: 1 for its size and 1 for its values
    l_ptr_size = sizeof(int) * (a_matrix->m + 1);

    memcpy((void *)l_free_address, &(l_ptr_size), sizeof(uint64_t));
    l_free_address += sizeof(uint64_t);

    // Align ptr adress on cache size
    l_padding = ( ( ( l_free_address - a_buffer_start_address + l_alignment - 1 ) / l_alignment ) * l_alignment ); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, l_csr->ptr, l_ptr_size);
    l_free_address += l_ptr_size;

    // Processing ind field.
    // 2 field are created: 1 for its size and 1 for its values
    l_ind_size = sizeof(int) * ( l_csr->ptr[a_matrix->m] - l_csr->base_index );

    // Align on 64bits since previous data are a multiple of int size
    l_padding = (( (l_free_address - a_buffer_start_address + 7) / 8 ) * 8); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, &(l_ind_size), sizeof(uint64_t));
    l_free_address += sizeof(uint64_t);

    // Align ind adress on cache size
    l_padding = (( (l_free_address - a_buffer_start_address + l_alignment - 1) / l_alignment ) * l_alignment); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, l_csr->ind, l_ind_size);
    l_free_address += l_ind_size;

    // Processing val field.
    // 2 field are created: 1 for its size and 1 for its values
    if ( a_matrix->type_value == COMPLEX_VALUE ) {
        l_val_size = sizeof(double) * 2 * ( l_csr->ptr[a_matrix->m] - l_csr->base_index );
    } else {
        l_val_size = sizeof(double) * ( l_csr->ptr[a_matrix->m] - l_csr->base_index );
    }

    // Align on 64bits since previous data are a multiple of int size
    l_padding = ( ( ( l_free_address - a_buffer_start_address + 7 ) / 8 ) * 8); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, &(l_val_size), sizeof(uint64_t));
    l_free_address += sizeof(uint64_t);

    // Align val adress on cache size
    l_padding = (( (l_free_address - a_buffer_start_address + l_alignment - 1) / l_alignment ) * l_alignment); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, l_csr->val, l_val_size);
    l_free_address += l_val_size;

    l_padding = (( (l_free_address - a_buffer_start_address + 7) / 8 ) * 8); 
    l_free_address = a_buffer_start_address + l_padding;

    memcpy((void *)l_free_address, &(l_csr->is_shared), sizeof(int));
    l_free_address += sizeof(uint64_t);

    // Update the address provided by caller
    *a_free_address = l_free_address;
}

/**
 * @brief 
 * 
 * @param a_address 
 * @return dmatCSR_t 
 */
dmatCSR_t VRP_Matrix_CSR_serializer::fromBuffer(uint64_t * a_address, uint64_t a_buffer_start_address) {
    dmatCSR_t l_csr  = (dmatCSR_t) malloc(sizeof(_dmatCSR_t));
    uint64_t l_next_address = *a_address;
    uint64_t l_ptr_size, l_ind_size, l_val_size;
    uint64_t l_alignment = VRP_Matrix_CSR_serializer::getAlignment();
    uint64_t l_padding = 0;

    if ( l_csr == NULL) {
        std::cout << "Fail allocating memory for oski_matCSR_t structure." << std::endl;
        return NULL;
    }

    l_csr->base_index = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_csr->has_unit_diag_implicit = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_csr->has_sorted_indices = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_csr->stored.is_upper = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_csr->stored.is_lower = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    // Processing ptr field.
    // 2 field are created: 1 for its size and 1 for its values
    l_ptr_size = *(uint64_t *)l_next_address;
    l_next_address += sizeof(uint64_t);

    // Align ptr adress on cache size
    l_padding = ( ( ( l_next_address - a_buffer_start_address + l_alignment - 1 ) / l_alignment ) * l_alignment ); 
    l_next_address = a_buffer_start_address + l_padding; 

    l_csr->ptr = (int *)l_next_address;
    l_next_address += l_ptr_size;

    l_padding = ( ( ( l_next_address - a_buffer_start_address + 7 ) / 8 ) * 8); 
    l_next_address = a_buffer_start_address + l_padding;

    // Processing ind field.
    // 2 field are created: 1 for its size and 1 for its values
    l_ind_size = *(uint64_t *)l_next_address;

    l_next_address += sizeof(uint64_t);

    // Align ind adress on cache size
    l_padding = ( ( ( l_next_address - a_buffer_start_address + l_alignment - 1 ) / l_alignment ) * l_alignment ); 
    l_next_address = a_buffer_start_address + l_padding; 

    l_csr->ind = (int *)l_next_address;
    l_next_address += l_ind_size;

    l_padding = ( ( ( l_next_address - a_buffer_start_address + 7 ) / 8 ) * 8 ); 
    l_next_address = a_buffer_start_address + l_padding;

    // Processing val field.
    // 2 field are created: 1 for its size and 1 for its values
    l_val_size = *(uint64_t *)l_next_address;
    l_next_address += sizeof(uint64_t);

    // Align val adress on cache size
    l_padding = ( ( ( l_next_address - a_buffer_start_address + l_alignment - 1 ) / l_alignment ) * l_alignment ); 
    l_next_address = a_buffer_start_address + l_padding; 

    l_csr->val = (double *)l_next_address;
    l_next_address += l_val_size;

    l_padding = ( ( ( l_next_address - a_buffer_start_address + 7 ) / 8 ) * 8); 
    l_next_address = a_buffer_start_address + l_padding;

    l_csr->is_shared = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    // Update the address provided by caller
    *a_address = l_next_address;

    return l_csr;
}

/**
 * @brief 
 * 
 * @param a_repr 
 */
void VRP_Matrix_CSR_serializer::print(matrix_t a_matrix) {
    dmatCSR_t l_csr = (dmatCSR_t)a_matrix->matrix->repr;
    if ( l_csr ==  NULL ) {
        printf("NULL\n");
        return;
    }

	std::cout << "base_index : "<< l_csr->base_index << std::endl;
    std::cout << "has_unit_diag_implicit : " << l_csr->has_unit_diag_implicit << std::endl;
    std::cout << "has_sorted_indices : " << l_csr->has_sorted_indices << std::endl;
    std::cout << "stored.is_upper : " << l_csr->stored.is_upper << std::endl;
    std::cout << "stored.is_lower : " <<  l_csr->stored.is_lower << std::endl;

    std::cout << "ptr : ";
    for (int l_ptr_index = 0 ; l_ptr_index < a_matrix->m + 1; l_ptr_index++) {
        std::cout << l_csr->ptr[l_ptr_index] << " ";
    }
    std::cout << std::endl;

    std::cout << "ind : ";
    for (int l_ind_index = 0 ; l_ind_index < l_csr->ptr[a_matrix->m] - l_csr->base_index; l_ind_index++) {
        std::cout << l_csr->ind[l_ind_index] << " ";
    }
    std::cout << std::endl;

    std::cout << "val : ";
    for (int l_val_index = 0 ; l_val_index < l_csr->ptr[a_matrix->m] - l_csr->base_index; l_val_index++) {
        std::cout << l_csr->val[l_val_index] << " - ";
    }
    std::cout << std::endl;

    std::cout << "is_shared : " << l_csr->is_shared << std::endl;
}
