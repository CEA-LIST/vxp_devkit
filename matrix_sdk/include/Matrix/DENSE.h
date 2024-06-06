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

#ifndef _DENSE_H_
#define _DENSE_H_

#include "complex.h"

/**
 *  OSKI DENSE format interface
 *  Double precision
 */
typedef struct __zmatDENSE_t {
    int lead_dim;
    complex_value_t * val;
} _zmatDENSE_t;

typedef _zmatDENSE_t * zmatDENSE_t;

/**
 *  OSKI DENSE format interface
 *  Double precision
 */
typedef struct __dmatDENSE_t {
    int lead_dim;
    double * val;
} _dmatDENSE_t;

typedef _dmatDENSE_t * dmatDENSE_t;

#endif /* _DENSE_H_ */
