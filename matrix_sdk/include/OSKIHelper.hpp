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

#ifndef __OSKIHELPER_HPP__
#define __OSKIHELPER_HPP__

#include "Matrix/matrix.h"

extern "C" {
    #include "oski/oski_Tid.h"
    // Reset bindings
    #define OSKI_REBIND 
    #include "oski/oski_Tiz.h"
}

typedef struct {
    bool complex;
    union {
        oski_matrix_t_Tid real_matrix;
        oski_matrix_t_Tiz complex_matrix;
    } oski_matrix;
} oski_matrix_wrapper_t;

namespace VPFloatPackage::OSKIHelper {  
    oski_matrix_wrapper_t loadFromFile(char * a_file_path);    

    oski_matrix_wrapper_t fromCSRMatrix(matrix_t a_matrix);

    matrix_t toBCSR(oski_matrix_wrapper_t a_sparse_matrix, int a_block_row_size, int a_block_col_size);

    oski_matrix_wrapper_t transpose(oski_matrix_wrapper_t a_matrix);

    // ATTENTION : OSKI convert DENSE matrices to Column Major format
    matrix_t toDense(oski_matrix_wrapper_t a_matrix, bool a_transpose=false,int a_lda=0);

    matrix_t toMatrix(oski_matrix_wrapper_t a_oski_matrix, char a_format=MATRIX_ROW_MAJOR, int a_lda=0, bool a_tuned_mat_flag=false);

    oski_matrix_wrapper_t buildCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, double * a_val, int a_base_index);

    void * getDiag(matrix_t a_input_matrix);

    extern bool oski_helper_initialized;
}

#endif /* __OSKIHELPER_HPP__ */