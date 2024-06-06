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

#include "VRPOffload/vrp_Matrix_BCSR_serializer.hpp"

using namespace VPFloatPackage::Offloading;

#define ALIGN_ADDRESS_64Bits(__offset) __offset = ( ( ( __offset / 8 ) + 1  ) * 8)

size_t VRP_Matrix_BCSR_serializer::getSizeBCSR(dmatBCSR_t a_bcsr, types_value_e a_type_value) {
    size_t l_size = 0;
    int l_total_nb_blocks = 0;

    if ( a_bcsr != NULL ) {
        // Add size for has_unit_diag_implicit field
        l_size += sizeof(uint64_t);

        // Add size for row_block_size, col_block_size, num_block_rows, num_block_cols fields
        l_size += sizeof(uint64_t) * 4;

        // Add size for bptr field
        l_size += sizeof(int) * (a_bcsr->num_block_rows + 1);
        if ( l_size / 8 ) {
            l_size += ( ( (l_size / 8 ) + 1 ) * 8 ) - l_size;
        }

        for ( int l_row_ptr = 0; l_row_ptr < a_bcsr->num_block_rows ; l_row_ptr++ ) {
            l_total_nb_blocks += a_bcsr->bptr[l_row_ptr+1] - a_bcsr->bptr[l_row_ptr];
        }

        // Add size for bind field
        l_size += sizeof(int) * l_total_nb_blocks;
        if ( l_size / 8 ) {
            l_size += ( ( (l_size / 8 ) + 1 ) * 8 ) - l_size;
        }

        // Add size for bval;
        if ( a_type_value == COMPLEX_VALUE ) {
            l_size += sizeof(double) * 2 * l_total_nb_blocks *  a_bcsr->row_block_size * a_bcsr->col_block_size;
        } else {
            l_size += sizeof(double) * l_total_nb_blocks *  a_bcsr->row_block_size * a_bcsr->col_block_size;
        }

        if ( l_size / 8 ) {
            l_size += ( ( (l_size / 8 ) + 1 ) * 8 ) - l_size;
        }

        // Add size for num_rows_leftover field
        l_size += sizeof(uint64_t);

        // Ignore size of mod_cached field since it is unused
        if ( a_bcsr->num_rows_leftover == 0 ) {
            l_size += sizeof(void *);
        } else {
            printf("Recursive size count engaged.\n");
            l_size += getSizeBCSR((dmatBCSR_t)(a_bcsr->leftover), a_type_value);
        }

        // Add size for mod_name field
        l_size += ( strlen(a_bcsr->mod_name) + 1 )  * sizeof(char);
        if ( l_size / 8 ) {
            l_size += ( ( (l_size / 8 ) + 1 ) * 8 ) - l_size;
        }
        // Add size for mod_cached
        l_size += sizeof(void *);
    }

    return l_size;
}


/**
 * @brief 
 * 
 * @param a_repr 
 * @return size_t 
 */
size_t VRP_Matrix_BCSR_serializer::getSize(matrix_t a_matrix) {
    dmatBCSR_t l_bcsr = (dmatBCSR_t)a_matrix->matrix->repr;
    return getSizeBCSR(l_bcsr, a_matrix->type_value);
}

/**
 * @brief 
 * 
 * @param a_bcsr 
 * @param a_address 
 */
void VRP_Matrix_BCSR_serializer::flatenBCSR(dmatBCSR_t a_bcsr, types_value_e a_type_value, uint64_t * a_free_address) {
    uint64_t l_free_address = *a_free_address;
    int l_total_nb_blocks = 0;
    uint64_t l_bval_size = 0;

    memcpy((void *)l_free_address, &(a_bcsr->has_unit_diag_implicit), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_bcsr->row_block_size), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_bcsr->col_block_size), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_bcsr->num_block_rows), sizeof(int));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, &(a_bcsr->num_block_cols), sizeof(int));
    l_free_address += sizeof(uint64_t);

    for ( int l_row_ptr = 0; l_row_ptr < a_bcsr->num_block_rows ; l_row_ptr++ ) {
        l_total_nb_blocks += a_bcsr->bptr[l_row_ptr+1] - a_bcsr->bptr[l_row_ptr];
    }

    memcpy((void *)l_free_address, a_bcsr->bptr, sizeof(int) * ( a_bcsr->num_block_rows + 1 ) );
    l_free_address += sizeof(int) * ( a_bcsr->num_block_rows + 1 );
    ALIGN_ADDRESS_64Bits(l_free_address);

    memcpy((void *)l_free_address, a_bcsr->bind, sizeof(int) * l_total_nb_blocks );
    l_free_address += sizeof(int) * l_total_nb_blocks;
    ALIGN_ADDRESS_64Bits(l_free_address);

    if ( a_type_value == COMPLEX_VALUE ) {
        l_bval_size = sizeof(double) * 2 * ( l_total_nb_blocks * a_bcsr->row_block_size * a_bcsr->col_block_size );
    } else {
        l_bval_size = sizeof(double) * ( l_total_nb_blocks * a_bcsr->row_block_size * a_bcsr->col_block_size );
    }
    
    memcpy((void *)l_free_address, &(l_bval_size), sizeof(uint64_t));
    l_free_address += sizeof(uint64_t);

    memcpy((void *)l_free_address, a_bcsr->bval, l_bval_size);
    l_free_address += l_bval_size;
    ALIGN_ADDRESS_64Bits(l_free_address);

    memcpy((void *)l_free_address, &(a_bcsr->num_rows_leftover), sizeof(int));
    l_free_address += sizeof(uint64_t);

    if ( a_bcsr->num_rows_leftover == 0 ) {
        l_free_address += sizeof(void *);
    } else {
        flatenBCSR(a_bcsr->leftover, a_type_value, &l_free_address);
    }

    strcpy((char *)l_free_address, a_bcsr->mod_name);
    l_free_address += sizeof(char) * strlen(a_bcsr->mod_name) + 1;
    ALIGN_ADDRESS_64Bits(l_free_address);

    l_free_address += sizeof(void *);

    // Update the address provided by caller
    *a_free_address = l_free_address;
}


/**
 * @brief 
 * 
 * @param a_matrix 
 * @param a_address 
 */
void VRP_Matrix_BCSR_serializer::flaten(matrix_t a_matrix, uint64_t * a_free_address) {
    dmatBCSR_t l_bcsr = (dmatBCSR_t)a_matrix->matrix->repr;
    uint64_t l_free_address = *a_free_address;

    flatenBCSR(l_bcsr, a_matrix->type_value, &l_free_address);

    // Update the address provided by caller
    *a_free_address = l_free_address;
}

/**
 * @brief 
 * 
 * @param a_address 
 * @return dmatBCSR_t 
 */
dmatBCSR_t VRP_Matrix_BCSR_serializer::fromBuffer(uint64_t * a_address) {
    dmatBCSR_t l_bcsr  = (dmatBCSR_t) malloc(sizeof(_dmatBCSR_t));
    uint64_t l_next_address = *a_address;
    int l_total_nb_blocks = 0;
    uint64_t l_bval_size = 0;

    if ( l_bcsr == NULL) {
        std::cout << "Fail allocating memory for oski_matBCSR_t structure." << std::endl;
        return NULL;
    }

    l_bcsr->has_unit_diag_implicit = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->row_block_size = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->col_block_size = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->num_block_rows = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->num_block_cols = *(int *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->bptr = (int*)l_next_address;
    l_next_address += sizeof(int) * ( l_bcsr->num_block_rows + 1 );
    ALIGN_ADDRESS_64Bits(l_next_address);

    for ( int l_row_ptr = 0; l_row_ptr < l_bcsr->num_block_rows ; l_row_ptr++ ) {
        l_total_nb_blocks += l_bcsr->bptr[l_row_ptr+1] - l_bcsr->bptr[l_row_ptr];
    }

    l_bcsr->bind = (int*)l_next_address;
    l_next_address += sizeof(int) * l_total_nb_blocks;
    ALIGN_ADDRESS_64Bits(l_next_address);

    l_bval_size = *(uint64_t *)l_next_address;
    l_next_address += sizeof(uint64_t);

    l_bcsr->bval = (double *)l_next_address;
    l_next_address += l_bval_size;
    ALIGN_ADDRESS_64Bits(l_next_address);

    l_bcsr->num_rows_leftover = *(int*)l_next_address;
    l_next_address += sizeof(uint64_t);

    if ( l_bcsr->num_rows_leftover == 0 ) {
        l_next_address += sizeof(void *);
    } else {
        l_bcsr->leftover = fromBuffer(&l_next_address);
    }

    l_bcsr->mod_name = (char *)l_next_address;
    l_next_address += sizeof(char) * strlen(l_bcsr->mod_name) + 1;
    ALIGN_ADDRESS_64Bits(l_next_address);

    l_next_address += sizeof(void *);

    // Update the address provided by caller
    *a_address = l_next_address;

    return l_bcsr;
}

void VRP_Matrix_BCSR_serializer::printBCSR(dmatBCSR_t a_bcsr) {
    if ( a_bcsr ==  NULL ) {
        printf("NULL\n");
        return;
    }

	std::cout << "row_block_size : "<< a_bcsr->row_block_size << std::endl;
    std::cout << "col_block_size : " << a_bcsr->col_block_size << std::endl;
    std::cout << "num_block_rows : " << a_bcsr->num_block_rows << std::endl;
    std::cout << "num_block_cols : " << a_bcsr->num_block_cols << std::endl;
    std::cout << "mod_name : " <<  a_bcsr->mod_name << std::endl;
    std::cout << "num_rows_leftover : " << a_bcsr->num_rows_leftover << std::endl;

    int l_bind_global_offset = 0;
    int l_bval_global_offset = 0;

    for(int l_row_id = 0; l_row_id < a_bcsr->num_block_rows; l_row_id++ ) {
        int l_nb_blocks_in_row = a_bcsr->bptr[l_row_id+1] - a_bcsr->bptr[l_row_id];

        for(int l_block_offset = 0; l_block_offset < l_nb_blocks_in_row-1; l_block_offset++) {

            std::cout << "Block(" << l_row_id << ", "<< a_bcsr->bind[l_bind_global_offset] << ")" << std::endl;

            for( int l_block_row_index = 0; l_block_row_index < a_bcsr->row_block_size; l_block_row_index++ ) {
               for( int l_block_col_index = 0 ; l_block_col_index < a_bcsr->col_block_size; l_block_col_index++ ) {
                    std::cout << a_bcsr->bval[l_bval_global_offset];
                    l_bval_global_offset++;
               }
               std::cout << std::endl;
            }
            std::cout << std::endl;
            l_bind_global_offset++;
        }
    }

    if ( a_bcsr->num_rows_leftover != 0 ) {
        printBCSR(a_bcsr->leftover);
    }
}

/**
 * @brief 
 * 
 * @param a_repr 
 */
void VRP_Matrix_BCSR_serializer::print(matrix_t a_matrix) {
    dmatBCSR_t l_bcsr = (dmatBCSR_t)a_matrix->matrix->repr;
    printBCSR(l_bcsr);
}
