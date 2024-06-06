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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdlib.h>
#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif

#include "complex.h"

/**
 *  Macro to get the tuned matrix pointer of an OSKI matrix
 *  Parameter:     oski_matrix_t  oski_matrix
 *  Return:        oski_mat_t     matrix
 */
#define FROM_OSKI(oski_matrix) (&(oski_matrix)->tuned_mat)

/**
 *  OSKI types ID
 *  (do not change the order)
 */
typedef enum types {
    NOT_A_FORMAT,
    CSR,
    CSC,
    BCSR,
    CB,
    GCSR,
    DENSE,
    MBCSR,
    VBR
} types_e;

typedef enum types_value {
    REAL_VALUE,
    COMPLEX_VALUE
} types_value_e;

/**
 *  OSKI matrix structure interface
 *  (oski_matspecific_t)
 */
typedef struct __oski_mat_t {
    size_t type_id;
    void * repr;
} _oski_mat_t;

typedef _oski_mat_t * oski_mat_t;

#define MATRIX_ROW_MAJOR 'N'
#define MATRIX_COL_MAJOR 'C'

// #define MATRIX_OPTIMIZE_LDA(n) ( ( ( ( ( n / 8 ) + 6 ) / 7 ) * 7 ) * 8 )
#define MATRIX_OPTIMIZE_LDA(n) ( ( (n + ( 8 * 5 ) - 1 ) / (8 * 5) ) * 8 * 5 )

typedef struct __matrix_t {
    int m;
    int n;
    int base_index;
    int lda;
    char format;
    types_e type_matrix;
    types_value_e type_value;
    oski_mat_t matrix;
} _matrix_t;

typedef _matrix_t * matrix_t;

matrix_t buildCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, double * a_val, int a_base_index);
matrix_t buildComplexCSR(int a_num_rows, int a_num_cols, int * a_row_ptr, int * a_col_ind, complex_value_t * a_val, int a_base_index);
matrix_t buildDiagCSR(int a_num_rows, double * a_val, int a_base_index);

void displayMatrix(matrix_t a_matrix, int a_display_precision);
int displayMatrixCharacteristics(matrix_t a_matrix);
void displayCSRStruct(matrix_t a_matrix, int a_display_precision);

void freeMatrix(matrix_t a_matrix);

double get(matrix_t a_matrix, int a_m, int a_n);

#ifdef __cplusplus
}
#endif

#endif /* __MATRIX_H__ */