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

void freeDenseMatrix(_dmatDENSE_t * a_dense_matrix) {
    free(a_dense_matrix->val);
}

void freeCSRMatrix(_dmatCSR_t * a_dense_matrix) {
    std::cout << __FUNCTION__ << "Need to be implemented!!!!!" << std::endl;
}

void freeBCSRMatrix(_dmatBCSR_t * a_dense_matrix) {
    std::cout << __FUNCTION__ << "Need to be implemented!!!!!" << std::endl;
}

void freeMatrix(matrix_t a_matrix) {
    switch(a_matrix->matrix->type_id) {
        case DENSE:
            freeDenseMatrix((_dmatDENSE_t *)(a_matrix->matrix->repr));
            break;
        case CSR:
            freeCSRMatrix((_dmatCSR_t *)(a_matrix->matrix->repr));
            break;
        case BCSR:
            freeBCSRMatrix((_dmatBCSR_t *)(a_matrix->matrix->repr));
            break;
        default :
            std::cout << "type_id" << a_matrix->matrix->type_id << " not supported for display." << std::endl;
            break;
    }
    free(a_matrix->matrix->repr);
    free(a_matrix);
}


matrix_t buildCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, double * a_val, int a_base_index) {
    matrix_t l_matrix = NULL;

    l_matrix = (matrix_t)malloc(sizeof(_matrix_t));

    if ( l_matrix == NULL ) {
        std::cout << "Fail allocating memory for new _matrix_t structure." << std::endl; 
        return NULL;
    }

    l_matrix->m = a_num_rows;
    l_matrix->n = a_num_cols;
    l_matrix->base_index = a_base_index;
    l_matrix->type_matrix = CSR;
    l_matrix->type_value = REAL_VALUE;
    l_matrix->format = MATRIX_ROW_MAJOR;
    l_matrix->lda = MATRIX_OPTIMIZE_LDA(l_matrix->n);

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
    dmatCSR_t l_csr_matrix = (dmatCSR_t)(l_matrix->matrix->repr);

    l_csr_matrix->base_index = a_base_index;
    l_csr_matrix->has_unit_diag_implicit = 0;
    l_csr_matrix->has_sorted_indices = 0;
    l_csr_matrix->has_sorted_indices = 0;
    l_csr_matrix->stored.is_upper = 1;
    l_csr_matrix->stored.is_upper = 1;
    l_csr_matrix->is_shared = 1;

    l_csr_matrix->ptr = a_row_ptr;
    l_csr_matrix->ind = a_col_ind;
    l_csr_matrix->val = a_val;

    return l_matrix;
}

matrix_t buildDiagCSR(int a_num_rows, double * a_val, int a_base_index) {
    int * l_ptr = (int *)malloc(sizeof(int) * a_num_rows);
    int * l_ind = (int *)malloc(sizeof(int) * a_num_rows);

    for ( int i = 0; i< a_num_rows; i++) {
        l_ptr[i] = i;
        l_ind[i] = i + a_base_index;
    }

    return buildCSR(a_num_rows, a_num_rows, l_ptr, l_ind, a_val, a_base_index);
}

void displayDENSEMatrix(matrix_t a_matrix) {

    dmatDENSE_t l_dense_matrix = (dmatDENSE_t)(a_matrix->matrix->repr);

    std::cout << "lead_dim:" << l_dense_matrix->lead_dim << std::endl;

    int l_row, l_col;
    for(l_row=0; l_row < a_matrix->m; l_row++) {
        std::cout << std::setw(2) << l_row + a_matrix->base_index << " : ";
        for(l_col=0; l_col < a_matrix->n; l_col++) {
            std::cout << std::setw(6) << l_dense_matrix->val[l_row  * a_matrix->lda + l_col] <<" " ;
        }
        std::cout << "\n";
    }
}

void displaydmatBCSR(dmatBCSR_t a_bcsr, int a_m, int a_n, int a_start_row_number, double * a_dense_matrix_display) {
    int l_nb_elements_per_block = a_bcsr->row_block_size * a_bcsr->col_block_size;
    int l_block_index = 0;
    int l_block_val_index = 0;

    // Processing of row of blocks
    for ( int l_row_block = 0 ; l_row_block < a_bcsr->num_block_rows; l_row_block++ ) {
        
        int l_real_start_row_index = l_row_block * a_bcsr->row_block_size + a_start_row_number;

        // Processing a block
        for (int l_block = a_bcsr->bptr[l_row_block]; l_block < a_bcsr->bptr[l_row_block+1]; l_block++) {

            int l_real_start_col_index = a_bcsr->bind[l_block_index];

            // Process each block line
            for ( int l_row_in_block = 0; l_row_in_block < a_bcsr->row_block_size; l_row_in_block++ ) {

                int l_real_row_index = l_real_start_row_index + l_row_in_block;

                // Process each column in a block line
                for ( int l_col_in_block = 0; l_col_in_block < a_bcsr->col_block_size; l_col_in_block++ ) {

                    int l_element_in_block_offset = l_row_in_block * a_bcsr->col_block_size + l_col_in_block;

                    int l_real_col_index = l_real_start_col_index + l_col_in_block;
                    
                    a_dense_matrix_display[l_real_row_index * a_n + l_real_col_index] = double( a_bcsr->bval[l_block_val_index + l_element_in_block_offset] );
                }

            }

            l_block_val_index += l_nb_elements_per_block;
            l_block_index++;

        }
    }

    if ( a_bcsr->num_rows_leftover > 0 ) {
        displaydmatBCSR(a_bcsr->leftover, 
                        a_m, 
                        a_n, 
                        a_start_row_number + ( a_bcsr->num_block_rows * a_bcsr->row_block_size),
                        a_dense_matrix_display
                    );
    }
}

void displayBCSRMatrix(matrix_t a_matrix) { 
    double * l_dense_matrix_display = (double *)malloc(sizeof(double) * a_matrix->m * a_matrix->n);

    if ( l_dense_matrix_display == NULL ){
        std::cout << "Fail allocating memory for BCSR matrix display." << std::endl;
        return;
    }

    memset(l_dense_matrix_display, 0.0, sizeof(double) * a_matrix->m * a_matrix->n);

    dmatBCSR_t l_bcsr_matrix = (dmatBCSR_t)(a_matrix->matrix->repr);

    displaydmatBCSR(l_bcsr_matrix, a_matrix->m, a_matrix->n, 0, l_dense_matrix_display);

    for(int l_row_index = 0 ; l_row_index < a_matrix->m; l_row_index++){
        std::cout << std::setw(2) << l_row_index  + a_matrix->base_index << " : ";

        for(int l_col_index = 0 ; l_col_index < a_matrix->m; l_col_index++){
            
                std::cout << std::setw(6) << *(l_dense_matrix_display + ( l_row_index * a_matrix->n + l_col_index) ) << " ";
        }

        std::cout << std::endl;
    }

    free(l_dense_matrix_display);
}

void displayCSRMatrix(matrix_t a_matrix,int a_display_precision) {
    int l_row_index, l_col, l_current_col_index;

    dmatCSR_t l_csr_matrix = (dmatCSR_t)(a_matrix->matrix->repr);

    l_current_col_index = 0;

    // Iterate over rows
    for(l_row_index=0; l_row_index < a_matrix->m; l_row_index++) {
        std::cout << std::setw(4) << l_row_index + l_csr_matrix->base_index << " : " ;

        int l_col_count_for_current_line =l_csr_matrix->ptr[l_row_index + 1] - l_csr_matrix->ptr[l_row_index];

        // Check if line is empty
        if ( l_col_count_for_current_line != 0 ) {
             
            for (l_col = l_csr_matrix->base_index; l_col < a_matrix->n + l_csr_matrix->base_index ; l_col++ ) {
                if ( l_csr_matrix->ind[l_current_col_index] == l_col && l_col_count_for_current_line != 0 ) {
                    if ( a_display_precision != 0 ) {
                        std::cout << std::setw(6) << l_csr_matrix->val[l_current_col_index]  << " - ";
                    } else {
                        std::cout << "+" ;
                    }
                    l_current_col_index++;
                    l_col_count_for_current_line--;
                } else {
                    if ( a_display_precision != 0 ) {
                        std::cout << std::setw(6) << 0.0  << " - " ;
                    } else {
                        std::cout << " ";
                    }
                }
            }
            
        } else {
            // Line is empty, display only 0
            for (l_col = 0; l_col < a_matrix->n; l_col++ ) {
                if ( a_display_precision != 0 ) {
                    std::cout << std::setw(6) << 0.0  << " - ";
                } else {
                    std::cout << " ";
                }
            }
        }

        std::cout << std::endl;
    }
}

void displayCSRStruct(matrix_t a_matrix, int a_display_precision) {
    dmatCSR_t l_csr_matrix = (dmatCSR_t)(a_matrix->matrix->repr);

    int * l_ptr = (int *)l_csr_matrix->ptr;
    int * l_ind = (int *)l_csr_matrix->ind;
    double * l_val = (double *)l_csr_matrix->val;

    int l_previous_ptr = 0;

    for ( int i = 0; i< a_matrix->m; i++) {
        std::cout << "PTR : " << l_ptr[i] << std::endl;
        for ( int j = l_previous_ptr; j < l_ptr[i]; j++ ) {
            std::cout << "IND : " << l_ind[j] << " - VAL : " << l_val[j] << std::endl;
        }
        l_previous_ptr=l_ptr[i];
    }

}

int displayMatrixCharacteristics(matrix_t a_matrix) {
    switch( a_matrix->type_matrix ) {
      case DENSE:
        std::cout << "dense";
        break;
      case CSR:
        std::cout << "sparse CSR";
        break;
      case BCSR:
        std::cout << "sparse BCSR";
        break;
      default:
        std::cout << "Matrix type " << a_matrix->type_matrix << " is not known" << std::endl;
        return 1;
        break;
    }

    std::cout << " ";

    switch( a_matrix->type_value ) {
      case COMPLEX_VALUE:
        std::cout << "complex";
        break;
      case REAL_VALUE:
        std::cout << "real";
        break;
      default:
        std::cout << "Value type " << a_matrix->type_value << " is not known" << std::endl;
        return 1;
        break;
    }

    return 0;

}


void displayMatrix(matrix_t a_matrix, int a_display_precision) {
    std::streamsize l_default_precicion = std::cout.precision();
    std::cout << std::setprecision(a_display_precision);

    switch(a_matrix->matrix->type_id) {
        case DENSE:
            displayDENSEMatrix(a_matrix);
            break;
        case CSR:
            displayCSRMatrix(a_matrix, a_display_precision);
            break;
        case BCSR:
            displayBCSRMatrix(a_matrix);
            break;
        default :
            std::cout << "type_id" << a_matrix->matrix->type_id << " not supported for display." << std::endl;
            break;
    }

    std::cout << std::setprecision(l_default_precicion) << std::defaultfloat ;
}

double getDENSE(matrix_t a_matrix, int a_m, int a_n) {
    dmatDENSE_t l_dense_matrix = (dmatDENSE_t)(a_matrix->matrix->repr);

    if ( a_m <= a_matrix->m && a_n <= a_matrix->n ) {
        return *(l_dense_matrix->val + ( ( ( a_m - a_matrix->base_index) * a_matrix->lda + ( a_n - a_matrix->base_index ) ) ) );
    } else {
        return NAN;
    }
}

double getdmatBCSR(dmatBCSR_t a_bcsr, int a_m, int a_n, int a_start_row_number) {
    // Check  if row is in leftover
    if ( ( a_start_row_number + a_bcsr->row_block_size * a_bcsr->num_block_rows ) < a_m ) {
        return getdmatBCSR(a_bcsr->leftover,  a_m, a_n, a_start_row_number + ( a_bcsr->num_block_rows * a_bcsr->row_block_size));
    }

    int l_block_size = a_bcsr->row_block_size * a_bcsr->col_block_size;
    int l_row_block_index = ( a_m - a_start_row_number ) / a_bcsr->row_block_size;
    int l_block_row_start_index = l_row_block_index * a_bcsr->row_block_size;

    // Iterate over row block to retrieve the one containing the value
    for ( int l_block_index = a_bcsr->bptr[l_row_block_index]; l_block_index < a_bcsr->bptr[l_row_block_index + 1]; l_block_index++) {
        int l_block_col_start_index = a_bcsr->bind[l_block_index];

        // Check if the current block contains the desired position
        if ( a_n >= l_block_col_start_index && a_n < l_block_col_start_index + a_bcsr->col_block_size ) {

            // Compute the position of the value in the current from its absolute position
            int l_in_block_value_row_index = a_m - a_start_row_number - l_block_row_start_index;
            int l_in_block_value_col_index = a_n - l_block_col_start_index;

            int l_in_block_value_offset = l_in_block_value_row_index * a_bcsr->col_block_size + l_in_block_value_col_index;

            // Return the value
            return a_bcsr->bval[l_block_index * l_block_size + l_in_block_value_offset];
        }
    }

    // If the value is not found return NAN
    return NAN;
}

double getBCSR(matrix_t a_matrix, int a_m, int a_n) {
    dmatBCSR_t l_bcsr_matrix = (dmatBCSR_t)(a_matrix->matrix->repr);

    if ( a_m <= a_matrix->m && a_n <= a_matrix->n ) {
        return getdmatBCSR(l_bcsr_matrix, a_m - a_matrix->base_index, a_n - a_matrix->base_index, 0);
    } else {
        return NAN;
    }
}

double getCSR(matrix_t a_matrix, int a_m, int a_n) {
    dmatCSR_t l_csr_matrix = (dmatCSR_t)(a_matrix->matrix->repr);

    int l_ind_offset = 0;
    int l_row_count = 0;

    if ( a_m > a_matrix->m || a_n > a_matrix->n ) {
        return NAN;
    }

    for(int l_row_index = 0; l_row_index <= a_m; l_row_index++) {
        l_row_count++;
        if ( l_row_index + a_matrix->base_index == a_m ) {
            for(int l_ind_index = l_ind_offset; l_ind_index <= l_csr_matrix->ptr[l_row_index]; l_ind_index++ ){

                if ( l_csr_matrix->ind[l_ind_index] == a_n ) {
                    return l_csr_matrix->val[l_ind_index];
                } else if (l_csr_matrix->ind[l_ind_index] > a_n ) {
                    break;
                }
            }
        }

        l_ind_offset = l_csr_matrix->ptr[l_row_index];
    }

    return NAN;
}

double get(matrix_t a_matrix, int a_m, int a_n) {
    double l_value = NAN;

    if ( a_m >= a_matrix->base_index && a_n >= a_matrix->base_index ) {

        switch(a_matrix->matrix->type_id) {
            case DENSE:
                l_value = getDENSE(a_matrix, a_m, a_n);
                break;
            case CSR:
                l_value = getCSR(a_matrix, a_m, a_n);
                break;
            case BCSR:
                l_value = getBCSR(a_matrix, a_m, a_n);
                break;
            default:
                std::cout << "type_id" << a_matrix->matrix->type_id << " not supported for display." << std::endl;
                break;
        }

    }

    return l_value;
}