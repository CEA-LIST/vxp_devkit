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
 *  @file   compiler.h
 *  @author Cesar Fuguet Tortolero
 */
#ifndef __COMPILER_H__
#define __COMPILER_H__

#define __NOINLINE__      __attribute__ ((noinline))
#define __INLINE__        __attribute__ ((inline))
#define __ALWAYS_INLINE__ __attribute__ ((always_inline))
#define __align(x)        __attribute__ ((aligned(x)))
#define __PACKED__        __attribute__ ((packed))
#define __likely(x)       __builtin_expect((x),1)
#define __unlikely(x)     __builtin_expect((x),0)

#endif /* __COMPILER_H__ */
