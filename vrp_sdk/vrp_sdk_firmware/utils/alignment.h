
/**
* Copyright 2022 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        alignment.c
 *  @author      Jerome Fereyre
 */

#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define ALIGN_ADDRESS_64Bits(__offset) __offset = ( ( ( __offset + 7 ) / 8 ) * 8 )
#define ALIGN_ADDRESS_64Bytes(__offset) __offset = ( ( ( __offset + 63 ) / 64 ) * 64 )

#ifdef __cplusplus
}
#endif

#endif /* __ALIGNMENT_H__ */
