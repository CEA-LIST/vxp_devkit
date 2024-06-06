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

#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>

#include "Matrix/matrix.h"
#include "Matrix/CSR.h"
#include "Matrix/DENSE.h"
#include "Matrix/BCSR.h"

matrix_t buildComplexCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, complex_value_t * a_val, int a_base_index) {
    matrix_t l_matrix = NULL;

    l_matrix = (matrix_t)malloc(sizeof(_matrix_t));

    if ( l_matrix == NULL ) {
        std::cout << "Fail allocating memory for new _matrix_t structure." << std::endl; 
        return NULL;
    }

    l_matrix->m = a_num_rows;
    l_matrix->n = a_num_cols;
    l_matrix->base_index = a_base_index;
    
    // Initialize matrix field
    l_matrix->matrix = (oski_mat_t)malloc(sizeof(_oski_mat_t));

    if ( l_matrix->matrix == NULL ) {
        std::cout << "Fail allocating memory for _oski_mat_t structure." << std::endl; 
        free(l_matrix);
        return NULL;
    }

    l_matrix->matrix->type_id = CSR;

    // Initialize matrix repr field for CSR matrix storage 
    l_matrix->matrix->repr = (dmatCSR_t)malloc(sizeof(_dmatCSR_t));

    if ( l_matrix->matrix->repr == NULL ) {
        std::cout << "Fail allocating memory for _dmatCSR_t structure." << std::endl;
        free(l_matrix->matrix);
        free(l_matrix);
        return NULL;
    }

    // Initialize CSR matrix data
    zmatCSR_t l_complex_csr_matrix = (zmatCSR_t)(l_matrix->matrix->repr);

    l_complex_csr_matrix->base_index = a_base_index;
    l_complex_csr_matrix->has_unit_diag_implicit = 0;
    l_complex_csr_matrix->has_sorted_indices = 0;
    l_complex_csr_matrix->has_sorted_indices = 0;
    l_complex_csr_matrix->stored.is_upper = 1;
    l_complex_csr_matrix->stored.is_upper = 1;
    l_complex_csr_matrix->is_shared = 1;

    l_complex_csr_matrix->ptr = a_row_ptr;
    l_complex_csr_matrix->ind = a_col_ind;
    l_complex_csr_matrix->val = a_val;

    return l_matrix;    
}
