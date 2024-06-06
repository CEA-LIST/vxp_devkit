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

#include "OSKIHelper.hpp"
#include "MTXUtil/crs.h"
#include "MTXUtil/MTXParser.hpp"
#include "MTXUtil/mmio.h"
#include "Matrix/CSR.h"
#include "Matrix/BCSR.h"
#include "Matrix/DENSE.h"
#include <iostream>
#include <string.h>
#include <stdio.h>

namespace VPFloatPackage::OSKIHelper {

};

using namespace VPFloatPackage;

bool oski_helper_initialized = false;

void oski_helper_initialize() {
    oski_Init();
    
    oski_helper_initialized = true;
}

void oski_helper_finalize() {
    oski_Close();

    oski_helper_initialized = false;
}

oski_matrix_wrapper_t oski_helper_create_real_CSR_matrix(crs_t * a_csr) {
    oski_matrix_wrapper_t l_matrix;

    oski_helper_initialize();

    l_matrix.complex = false;
    l_matrix.oski_matrix.real_matrix = oski_CreateMatCSR_Tid (
        a_csr->ROW_PTR,
        a_csr->COL_IND,
        a_csr->VAL,
        a_csr->M,
        a_csr->N,
        SHARE_INPUTMAT,
        1,
        INDEX_ONE_BASED );

    oski_helper_finalize();

    return l_matrix;
}

oski_matrix_wrapper_t oski_helper_create_complex_CSR_matrix(crs_t * a_csr) {
    oski_matrix_wrapper_t l_matrix;

    oski_helper_initialize();

    l_matrix.complex = true;
    l_matrix.oski_matrix.complex_matrix = oski_CreateMatCSR_Tiz (
        a_csr->ROW_PTR,
        a_csr->COL_IND,
        (doublecomplex_t*)(a_csr->VAL),
        a_csr->M,
        a_csr->N,
        SHARE_INPUTMAT,
        1,
        INDEX_ONE_BASED );

    oski_helper_finalize();

    return l_matrix;
}

oski_matrix_wrapper_t OSKIHelper::loadFromFile(char * a_file_path) {
    crs_t * l_csr = ::MTXParser::parseFileToCRS(a_file_path);

    if ( l_csr == NULL ) {
        std::cout << __FUNCTION__ << "MTX file parsing failed. return NULL" << std::endl;
        oski_matrix_wrapper_t l_wrapper;
        l_wrapper.complex = false;
        l_wrapper.oski_matrix.real_matrix = NULL;
        return l_wrapper;
    }

    if ( l_csr->VALUE_TYPE == COMPLEX ) {
        return oski_helper_create_complex_CSR_matrix(l_csr);
    } else {
        return oski_helper_create_real_CSR_matrix(l_csr);
    }
}

types_e getOSKIMatTypeId(oski_id_t a_type, oski_id_t a_oski_ind_id, oski_id_t a_oski_val_id) {
    if ( a_type == oski_LookupMatTypeId("CSR", a_oski_ind_id, a_oski_val_id) ) {
        return CSR;
    }

    if ( a_type == oski_LookupMatTypeId("BCSR", a_oski_ind_id, a_oski_val_id) ) {
        return BCSR;
    }

    if ( a_type == oski_LookupMatTypeId("DENSE", a_oski_ind_id, a_oski_val_id) ) {
        return DENSE;
    }

    return NOT_A_FORMAT;
}

/**
 * Compute optimzed LDA value depending on number of element in each matrix line and their size
 * n                   : number of element in matrix line
 * matrix_element_size : should be power of 2, if not the case, we should rework
 */
int ComputeOptimizedLDA(int n, int matrix_element_size) {
    // size of cache line (assumed to be the same for all architectures)
    int cache_line_size = 64;

    // how many elements per cache line
    int elem_per_cache_line = cache_line_size / matrix_element_size; // integer as everything is power of 2

    // first align matrix LDA to multiple of cache line
    int LDA = ( ( ( n + elem_per_cache_line - 1 ) / elem_per_cache_line ) * elem_per_cache_line );

    // now compute the number of cache lines for given LDA:
    int LDA_cache_lines = ( LDA * matrix_element_size ) / cache_line_size; // should be integer

    // if LDA cache_lines is even, add 1 to get odd cache lines
    if( ! ( LDA_cache_lines & 1 ) ) {
        LDA_cache_lines += 1;

        // compute new LDA
        LDA = (LDA_cache_lines * cache_line_size) / matrix_element_size;
    }

    return LDA;
}

matrix_t OSKIHelper::toMatrix(oski_matrix_wrapper_t a_oski_matrix, char a_format, int a_lda, bool a_tuned_mat_flag) {
    matrix_t l_matrix =(matrix_t)malloc(sizeof(_matrix_t));

    if (l_matrix == NULL ) {
        std::cout << "Fail allocating memory for new _matrix_t strucuture." << std::endl;
        return NULL;
    }

    if ( a_oski_matrix.complex ) {
        l_matrix->m = a_oski_matrix.oski_matrix.complex_matrix->props.num_rows;
        l_matrix->n = a_oski_matrix.oski_matrix.complex_matrix->props.num_cols;
        l_matrix->base_index = 1;
        l_matrix->type_value = COMPLEX_VALUE;
        l_matrix->type_matrix = getOSKIMatTypeId(a_oski_matrix.oski_matrix.complex_matrix->tuned_mat.type_id, 1, 4);
        
        if ( a_tuned_mat_flag ) {
            l_matrix->matrix = (oski_mat_t)&(a_oski_matrix.oski_matrix.complex_matrix->tuned_mat);
        } else {
            l_matrix->matrix = (oski_mat_t)&(a_oski_matrix.oski_matrix.complex_matrix->input_mat);
        }
    } else {
        l_matrix->m = a_oski_matrix.oski_matrix.real_matrix->props.num_rows;
        l_matrix->n = a_oski_matrix.oski_matrix.real_matrix->props.num_cols;
        l_matrix->base_index = 1;
        l_matrix->type_value = REAL_VALUE;
        l_matrix->type_matrix = getOSKIMatTypeId(a_oski_matrix.oski_matrix.complex_matrix->tuned_mat.type_id, 1, 2);

        if ( a_tuned_mat_flag ) {
            l_matrix->matrix = (oski_mat_t)&(a_oski_matrix.oski_matrix.real_matrix->tuned_mat);    
        } else {
            l_matrix->matrix = (oski_mat_t)&(a_oski_matrix.oski_matrix.real_matrix->input_mat);
        }
    }

    l_matrix->format = a_format;

    if ( a_lda > 0 ) {
        l_matrix->lda = a_lda;
    } else {
        if ( l_matrix->format == 'N' ) {
            l_matrix->lda = ComputeOptimizedLDA(l_matrix->n, l_matrix->type_value == REAL_VALUE ? 8 : 16);
        } else {
            l_matrix->lda = ComputeOptimizedLDA(l_matrix->m, l_matrix->type_value == REAL_VALUE ? 8 : 16);
        }
    }

    return l_matrix;
}

oski_matrix_wrapper_t OSKIHelper::fromCSRMatrix(matrix_t a_matrix) {
    oski_matrix_wrapper_t l_oski_matrix;

    if ( a_matrix->matrix->type_id != CSR ) {
        std::cout << "Input matrix is not a CSR matrix." << std::endl;
        l_oski_matrix.oski_matrix.real_matrix = NULL;
        return l_oski_matrix;
    }

    dmatCSR_t l_csr_matrix = (dmatCSR_t)(a_matrix->matrix->repr);

    oski_helper_initialize();

    l_oski_matrix.oski_matrix.real_matrix = oski_CreateMatCSR_Tid (
        l_csr_matrix->ptr,
        l_csr_matrix->ind,
        l_csr_matrix->val,
        a_matrix->m,
        a_matrix->n,
        SHARE_INPUTMAT,
        1,
        INDEX_ONE_BASED );

    oski_helper_finalize();

    l_oski_matrix.complex = false;

    return l_oski_matrix;
}
matrix_t OSKIHelper::toDense(oski_matrix_wrapper_t a_sparse_matrix, bool a_transpose, int a_lda) {
    dmatCSR_t l_input_csr = NULL;
    uint8_t l_scalar_type_offset_factor=1;
    matrix_t l_dense_matrix = (matrix_t)malloc(sizeof(_matrix_t));
    l_dense_matrix->matrix = (oski_mat_t)malloc(sizeof(_oski_mat_t));    

    dmatDENSE_t l_matrix_data_struct = (dmatDENSE_t)malloc(sizeof(_dmatDENSE_t));

    l_dense_matrix->base_index = 1;
    l_dense_matrix->type_matrix = DENSE;
    l_dense_matrix->matrix->type_id = l_dense_matrix->type_matrix;

    // Deal with type of data stored within the matrix (real or complex saclars)
    if ( a_sparse_matrix.complex ) {
        l_dense_matrix->type_value = COMPLEX_VALUE;
        l_input_csr = (dmatCSR_t)(a_sparse_matrix.oski_matrix.complex_matrix->input_mat.repr);
        l_scalar_type_offset_factor=2;
    } else {
        l_dense_matrix->type_value = REAL_VALUE;
        l_input_csr = (dmatCSR_t)(a_sparse_matrix.oski_matrix.real_matrix->input_mat.repr);
        l_scalar_type_offset_factor=1;
    }

    if ( a_transpose ) {
        l_dense_matrix->n = a_sparse_matrix.oski_matrix.real_matrix->props.num_rows;
        l_dense_matrix->m = a_sparse_matrix.oski_matrix.real_matrix->props.num_cols;     
        
    } else {
        l_dense_matrix->m = a_sparse_matrix.oski_matrix.real_matrix->props.num_rows;
        l_dense_matrix->n = a_sparse_matrix.oski_matrix.real_matrix->props.num_cols;        
    }

    l_dense_matrix->format = MATRIX_ROW_MAJOR;
    if ( a_lda != 0 ) {
        l_dense_matrix->lda = a_lda;
        std::cout << "LDA set manually to " << l_dense_matrix->lda << std::endl;
    } else {
        // Compute LDA
        if ( l_dense_matrix->n == 1 ) {
            l_dense_matrix->lda = l_dense_matrix->n;
        } else {
            l_dense_matrix->lda = ComputeOptimizedLDA(l_dense_matrix->n, l_dense_matrix->type_value == REAL_VALUE ? 8 : 16);
        }
        std::cout << "LDA set automatically to " << l_dense_matrix->lda << std::endl;
    }

    l_matrix_data_struct->lead_dim = l_dense_matrix->n;

    // Allocate memory for val array
    size_t l_val_size = sizeof(double) * l_dense_matrix->m * l_dense_matrix->lda;
    l_matrix_data_struct->val = (double *)malloc(l_val_size);

    // Initialize array to 0
    memset(l_matrix_data_struct->val, 0, l_val_size);

    // Iterate over matrix row indeces.
    for (size_t l_row_index = 0; l_row_index < a_sparse_matrix.oski_matrix.real_matrix->props.num_rows; l_row_index++)
    {
        // Iterate over col indices for the current CSR row
        for (size_t l_csr_col_index = l_input_csr->ptr[l_row_index] - l_input_csr->base_index; l_csr_col_index < l_input_csr->ptr[l_row_index + 1] - l_input_csr->base_index; l_csr_col_index++)
        {
            size_t l_real_col_index = l_input_csr->ind[l_csr_col_index];
            size_t l_val_dense_val_offset = 0;

            // Compute the offset of the current value in the destination dense matrix.
            // If a transposition is required twek the coordinate accordingly.
            if ( a_transpose == false)  {
                l_val_dense_val_offset = (l_real_col_index - l_input_csr->base_index) + (l_row_index) * l_dense_matrix->lda * l_scalar_type_offset_factor;
            } else {
                l_val_dense_val_offset = (l_row_index) + (l_real_col_index- l_input_csr->base_index) * l_dense_matrix->lda * l_scalar_type_offset_factor;
            }

            // Store the current CSR value in the destination dense matrix.
            l_matrix_data_struct->val[l_val_dense_val_offset] += l_input_csr->val[l_csr_col_index ];

            if ( a_sparse_matrix.complex ) {
                l_matrix_data_struct->val[l_val_dense_val_offset+1] += l_input_csr->val[l_csr_col_index - l_dense_matrix->base_index + 1];
            }
        }
    }

    l_dense_matrix->matrix->repr = l_matrix_data_struct;

    return l_dense_matrix;
}

matrix_t OSKIHelper::toBCSR(oski_matrix_wrapper_t a_sparse_matrix, int a_block_row_size, int a_block_col_size) {
    char l_lua_transform[1024];

    matrix_t l_bcsr_matrix = (matrix_t)malloc(sizeof(_matrix_t));
    l_bcsr_matrix->base_index = 1;
    l_bcsr_matrix->matrix = (oski_mat_t)malloc(sizeof(_oski_mat_t));
    l_bcsr_matrix->matrix->type_id = BCSR;
    l_bcsr_matrix->type_matrix = BCSR;
    l_bcsr_matrix->m = a_sparse_matrix.oski_matrix.complex_matrix->props.num_rows;
    l_bcsr_matrix->n = a_sparse_matrix.oski_matrix.complex_matrix->props.num_cols; 

    snprintf(l_lua_transform, 1024, "A_new = BCSR(InputMat, %d, %d)\n return A_new\n", a_block_row_size, a_block_col_size);

    oski_helper_initialize();

    if ( a_sparse_matrix.complex ) {
        oski_matrix_t_Tiz l_oski_bcsr_matrix = oski_CopyMat_Tiz (a_sparse_matrix.oski_matrix.complex_matrix);    

        oski_ApplyMatTransforms_Tiz(l_oski_bcsr_matrix, l_lua_transform);

        l_bcsr_matrix->type_value = COMPLEX_VALUE;
        l_bcsr_matrix->matrix->repr = l_oski_bcsr_matrix->tuned_mat.repr;
    } else {
        oski_matrix_t_Tid l_oski_bcsr_matrix = oski_CopyMat_Tid (a_sparse_matrix.oski_matrix.real_matrix);    

        oski_ApplyMatTransforms_Tid(l_oski_bcsr_matrix, l_lua_transform);

        l_bcsr_matrix->type_value = REAL_VALUE;
        l_bcsr_matrix->matrix->repr = l_oski_bcsr_matrix->tuned_mat.repr;
    }

    l_bcsr_matrix->format = MATRIX_ROW_MAJOR;
    l_bcsr_matrix->lda = ComputeOptimizedLDA(l_bcsr_matrix->n, l_bcsr_matrix->type_value == REAL_VALUE ? 8 : 16);

    oski_helper_finalize();

    return l_bcsr_matrix;
}

void transpose(oski_matrix_wrapper_t a_input_matrix , oski_matrix_wrapper_t a_transposed_matrix){
    int l_input_M,  l_input_N, l_input_NZ;
    int l_transposed_M, l_transposed_N, l_transposed_NZ;
    dmatCSR_t l_input_csr = NULL;
    dmatCSR_t l_transposed_csr = NULL;

    a_transposed_matrix.complex = a_input_matrix.complex;

    if ( a_transposed_matrix.complex ) {
        l_input_csr = (dmatCSR_t)(a_input_matrix.oski_matrix.complex_matrix->input_mat.repr);
        l_transposed_csr = (dmatCSR_t)(a_transposed_matrix.oski_matrix.complex_matrix->input_mat.repr);

        l_transposed_M = a_input_matrix.oski_matrix.complex_matrix->props.num_rows;
        l_transposed_N = a_input_matrix.oski_matrix.complex_matrix->props.num_cols;
        l_transposed_NZ = a_input_matrix.oski_matrix.complex_matrix->props.num_nonzeros;

        a_transposed_matrix.oski_matrix.complex_matrix->input_mat.type_id = a_input_matrix.oski_matrix.complex_matrix->input_mat.type_id;
        a_transposed_matrix.oski_matrix.complex_matrix->props.num_rows = l_transposed_M;
        a_transposed_matrix.oski_matrix.complex_matrix->props.num_cols = l_transposed_N;
        a_transposed_matrix.oski_matrix.complex_matrix->props.num_nonzeros = l_transposed_NZ;
    } else {
        l_input_csr = (dmatCSR_t)(a_input_matrix.oski_matrix.real_matrix->input_mat.repr);
        l_transposed_csr = (dmatCSR_t)(a_transposed_matrix.oski_matrix.real_matrix->input_mat.repr);

        l_transposed_M = a_input_matrix.oski_matrix.real_matrix->props.num_rows;
        l_transposed_N = a_input_matrix.oski_matrix.real_matrix->props.num_cols;
        l_transposed_NZ = a_input_matrix.oski_matrix.real_matrix->props.num_nonzeros;

        a_transposed_matrix.oski_matrix.real_matrix->input_mat.type_id = a_input_matrix.oski_matrix.real_matrix->input_mat.type_id;
        a_transposed_matrix.oski_matrix.real_matrix->props.num_rows = l_transposed_M;
        a_transposed_matrix.oski_matrix.real_matrix->props.num_cols = l_transposed_N;
        a_transposed_matrix.oski_matrix.real_matrix->props.num_nonzeros = l_transposed_NZ;
    }
    
    free(l_transposed_csr->ind);
    free(l_transposed_csr->ptr);
    free(l_transposed_csr->val);

    l_transposed_csr->ptr = (int *)malloc(sizeof(int) * (l_transposed_M + 1));
    memset(l_transposed_csr->ptr, 0, sizeof(int) * (l_transposed_M + 1));

    l_transposed_csr->ind = (int *)malloc(sizeof(int) * l_transposed_NZ);
    memset(l_transposed_csr->ind, 0, sizeof(int) * l_transposed_NZ);

    /*
     * Val field allocation depends on the type of value stored in the matrix:
     *  - real : sizeof(double) * NZNUM
     *  - complex : sizeof(double) * NZNUM * 2
     */
    size_t l_transposed_t_val_size  = sizeof(double) * l_transposed_NZ;

    if ( a_transposed_matrix.complex ) {
        l_transposed_t_val_size *= 2;
    }

    l_transposed_csr->val = (double *)malloc(l_transposed_t_val_size);
    memset(l_transposed_csr->val, 0, l_transposed_t_val_size);

    /*
     * Create an array used to store the number of non null elements in each column 
     * of the transposed matrix row.
     * This array will be reuse later to compute the col_ind value of the next element
     * to store in a row.
     * This array contains M elements, one per transposed matrix row
     */
    int * l_transposed_row_col_ind_info = (int *)malloc(sizeof(int) * l_transposed_M);
    memset(l_transposed_row_col_ind_info, 0, sizeof(int) *l_transposed_M);

    /*
     * Walk through the input matrix in order to compute the number of element in each column.
     */
    for ( int l_row_index = 0 ; l_row_index < l_transposed_M; l_row_index++) {
        for ( int l_col_index = l_input_csr->ptr[l_row_index]; l_col_index < l_input_csr->ptr[l_row_index+1]; l_col_index++ ) {
            l_transposed_row_col_ind_info[l_input_csr->ind[l_col_index-1]-1]++;
        }
    }

    /*
     * Based on the previous computation loop we are able to create ROW_PTR of the transposed matrix.
     * Each l_transposed_row_col_ind_info element is the number of element a transposed row.
     * Then the l_transposed_row_col_ind_info is reused to store the col_ind of the next element to store in a row
     */
    l_transposed_csr->ptr[0] = 1;
    for ( int i = 1 ; i <= l_transposed_M; i++) {
        l_transposed_csr->ptr[i] = l_transposed_csr->ptr[i-1] + l_transposed_row_col_ind_info[i-1];
        // Store in l_transposed_row_col_ind_info the start index of the col_ind for the previous transposed row
        l_transposed_row_col_ind_info[i-1] = l_transposed_csr->ptr[i-1];
    }
    // Store in l_transposed_row_col_ind_info the start index of the col_ind for the last transposed row
    l_transposed_row_col_ind_info[l_transposed_M-1] = l_transposed_csr->ptr[l_transposed_M-1];
    

    /*
     * Walk through the intput matrix and 
     */
    for ( int l_row_index = 0 ; l_row_index < l_transposed_M; l_row_index++) {
        for ( int l_col_index = l_input_csr->ptr[l_row_index]; l_col_index < l_input_csr->ptr[l_row_index+1]; l_col_index++ ) {

            // build the position of the current element in the transposed matrix
            int l_new_col_index = l_row_index + 1;
            int l_new_row_index = l_input_csr->ind[l_col_index-1];

            // Retrieve the col_ind of the current element in the transposed matrix
            int l_transposed_col_ind = l_transposed_row_col_ind_info[l_new_row_index-1]-1;

            l_transposed_csr->ind[l_transposed_col_ind] = l_new_col_index;
            if ( a_transposed_matrix.complex ) {
                l_transposed_csr->val[l_transposed_col_ind * 2] = l_input_csr->val[(l_col_index-1) * 2 ];
                l_transposed_csr->val[l_transposed_col_ind * 2 + 1] = l_input_csr->val[(l_col_index-1) * 2 + 1];
            } else {
                l_transposed_csr->val[l_transposed_col_ind] = l_input_csr->val[l_col_index-1];
            }

            // Move the col_ind to the next element of the row in the transposed matrix
            l_transposed_row_col_ind_info[l_new_row_index-1]++;
        }
    }

    free(l_transposed_row_col_ind_info);

}

oski_matrix_wrapper_t OSKIHelper::transpose(oski_matrix_wrapper_t a_matrix) {
    oski_matrix_wrapper_t l_transposed_matrix;

    l_transposed_matrix.complex = false;

    if ( a_matrix.complex ) {
        l_transposed_matrix.complex = true;

        oski_helper_initialize();

        l_transposed_matrix.oski_matrix.complex_matrix = oski_CopyMat_Tiz (a_matrix.oski_matrix.complex_matrix);

        transpose(a_matrix, l_transposed_matrix);

        oski_helper_finalize();

    } else {
        oski_helper_initialize();

        l_transposed_matrix.oski_matrix.real_matrix = oski_CopyMat_Tid (a_matrix.oski_matrix.real_matrix);

        transpose(a_matrix, l_transposed_matrix);

        oski_helper_finalize();
    }

    return l_transposed_matrix;    
}

oski_matrix_wrapper_t OSKIHelper::buildCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, double * a_val, int a_base_index) {
	oski_matrix_wrapper_t l_matrix;

    oski_helper_initialize();

    if ( a_base_index == 1 ) {
	    l_matrix.oski_matrix.real_matrix = oski_CreateMatCSR_Tid( 	a_row_ptr, a_col_ind, a_val, a_num_rows, a_num_cols, SHARE_INPUTMAT, 1, INDEX_ONE_BASED );
    } else {
        l_matrix.oski_matrix.real_matrix = oski_CreateMatCSR_Tid( 	a_row_ptr, a_col_ind, a_val, a_num_rows, a_num_cols, SHARE_INPUTMAT, 1, INDEX_ZERO_BASED );
    }

    l_matrix.complex = false;

    oski_helper_finalize();

    return l_matrix;
}

void * OSKIHelper::getDiag(matrix_t a_input_matrix) {
    oski_matrix_wrapper_t l_matrix;
    void * l_diag = NULL;

    /* recuperer les valeurs diag */
    oski_matrix_wrapper_t l_input_matrix = fromCSRMatrix(a_input_matrix);
    
    if ( a_input_matrix->type_value == COMPLEX_VALUE ) {
        l_diag = ( doublecomplex_t * ) malloc(a_input_matrix->n * sizeof(complex_t));    
    } else {
        l_diag = ( double * ) malloc(a_input_matrix->n * sizeof(double));    
    }

    oski_Init(); 

    for (int jj=1; jj <= a_input_matrix->n; jj++) {

        if ( a_input_matrix->type_value == COMPLEX_VALUE ) {
            ( ( doublecomplex_t * ) l_diag )[jj-1] = oski_GetMatEntry_Tiz(l_input_matrix.oski_matrix.complex_matrix, jj, jj);
        } else {
            ( ( double * ) l_diag )[jj-1] = oski_GetMatEntry_Tid(l_input_matrix.oski_matrix.real_matrix, jj, jj);
        }
    }

    oski_Close();
    
    return l_diag;
}