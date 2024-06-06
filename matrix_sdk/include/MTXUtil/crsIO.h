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
 * Authors       : Yves Durand, Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#ifndef __YDCRSIO__
#define __YDCRSIO__

#ifdef __cplusplus
extern "C" {
#endif

#include "MTXUtil/crs.h"

int printVectorFile(char * filename, int N, double * VEKTOR ) ;
int printFFVectorFile(char * filename, int N, double * VEKTOR ) ;

int readVectorFile(char * filename, int N, double * VEKTOR ) ;

int readMMfile(char * filename, struct crs * SPM );

#ifdef __cplusplus
}
#endif


#endif
