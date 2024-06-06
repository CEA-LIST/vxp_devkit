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
 * Authors       : Valentin Isaac--Chassande
 * Creation Date : August, 2023
 * Description   : 
 **/

#ifndef _BCSR_H_
#define _BCSR_H_

#include "complex.h"

/**
 *  OSKI BCSR format interface
 *  Double precision
 */
typedef struct __zmatBCSR_t {
    int has_unit_diag_implicit; // Unused for now

    int row_block_size; // r
    int col_block_size; // c

    int num_block_rows;
    int num_block_cols; // Including the last column which can be filled with zeros

    int * bptr;
    int * bind;
    complex_value_t * bval;

    int num_rows_leftover;

    struct __dmatBCSR_t * leftover; // NULL if no leftover row

    char * mod_name; // Unused
    void * mod_cached; // Unused
} _zmatBCSR_t;

typedef _zmatBCSR_t * zmatBCSR_t;

/**
 *  OSKI BCSR format interface
 *  Double precision
 */
typedef struct __dmatBCSR_t {
    int has_unit_diag_implicit; // Unused for now

    int row_block_size; // r
    int col_block_size; // c

    int num_block_rows;
    int num_block_cols; // Including the last column which can be filled with zeros

    int * bptr;
    int * bind;
    double * bval;

    int num_rows_leftover;

    struct __dmatBCSR_t * leftover; // NULL if no leftover row

    char * mod_name; // Unused
    void * mod_cached; // Unused
} _dmatBCSR_t;

typedef _dmatBCSR_t * dmatBCSR_t;

#endif /* _BCSR_H_ */
