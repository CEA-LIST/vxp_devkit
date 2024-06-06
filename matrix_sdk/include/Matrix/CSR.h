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

#ifndef _CSR_H_
#define _CSR_H_

#include "complex.h"

/**
 *  OSKI CSR format interface
 *  Double precision
 */
typedef struct __zmatCSR_t {
    int base_index;

    int has_unit_diag_implicit;

    int has_sorted_indices;

    struct
    {
        int is_upper;
        int is_lower;
    } stored;

    int *ptr;
    int *ind;
    complex_value_t *val;
    int is_shared;
} _zmatCSR_t;

typedef _zmatCSR_t * zmatCSR_t;


/**
 *  OSKI CSR format interface
 *  Double precision
 */
typedef struct __dmatCSR_t {
    int base_index;

    int has_unit_diag_implicit;

    int has_sorted_indices;

    struct
    {
        int is_upper;
        int is_lower;
    } stored;

    int *ptr;
    int *ind;
    double *val;
    int is_shared;
} _dmatCSR_t;

typedef _dmatCSR_t * dmatCSR_t;

#endif /* _CSR_H_ */
