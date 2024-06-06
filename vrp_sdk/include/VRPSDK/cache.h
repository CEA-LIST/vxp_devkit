/**
* Copyright 2021 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        cache.h
 *  @author      Cesar Fuguet-Tortolero
 *  @date        October, 2021
 */
#ifndef __CACHE_H__
#define __CACHE_H__

#ifdef __cplusplus
extern "C" {
#endif

#define VRP_CACHE_LINE_SIZE_IN_BYTES 64 

#define ALIGN_ADDRESS_64Bytes(__address) ( ( ( __address + VRP_CACHE_LINE_SIZE_IN_BYTES - 1 ) / VRP_CACHE_LINE_SIZE_IN_BYTES ) * VRP_CACHE_LINE_SIZE_IN_BYTES )

#ifdef __cplusplus
}
#endif

#endif /* __CACHE_H__ */
