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
 *  @file        cpu.c
 *  @author      Jerome Fereyre
 */
#include "VRPSDK/perfcounters/cpu.h"

#include <stdint.h>

#if defined(__x86_64__) || defined(_M_X64)
#include <x86intrin.h>
#endif /* __x86_64__ */


uint64_t get_cycles() {

#if defined(__x86_64__) || defined(_M_X64)
    /**********************************************************************************
     *                                      X86_64 
     *********************************************************************************/
    return __rdtsc();

#elif defined(__riscv)
    /**********************************************************************************
     *                                      RISCV 
     *********************************************************************************/
#if (__SIZEOF_LONG__ == 4)
	register uint32_t lo;
	register uint32_t hi;

	asm volatile (
			"csrr %[lo], mcycle                  \n"
			"csrr %[hi], mcycleh                 \n"
			: [lo] "=r"(lo), [hi] "=r"(hi)
			:
			: "memory");

	/* as there is one cycle between the reading of the lower bits and
	 * higher bits, decrement hi if lo was about to wrap */
	if (lo == 0xffffffffUL)
		hi--;

	return ((uint64_t)hi << 32) | lo;
#else
	register uint64_t cycles;

	asm volatile (
			"csrr %[cycles], mcycle            \n"
			: [cycles] "=r"(cycles)
			:
			: "memory");

	return cycles;
#endif /* __SIZEOF_LONG__ */

#else
    /**********************************************************************************
     *                                      UNSUPPORTED
     *********************************************************************************/
    return 0;

#endif
}

uint64_t cpu_imiss()
{
#if defined(__x86_64__) || defined(_M_X64)
	return 0;

#elif defined(__riscv)
#if (__SIZEOF_LONG__ == 4)
    register uint32_t lo;
    register uint32_t hi;

    asm volatile (
            "csrr %[lo], mhpmcounter3       \n"
            "csrr %[hi], mhpmcounter3h      \n"
            : [lo] "=r"(lo), [hi] "=r"(hi)
            :
            : "memory");

    /* as there is one cycle between the reading of the lower bits and
     * higher bits, decrement hi if lo was about to wrap */
    if (lo == 0xffffffffUL)
        hi--;

    return ((uint64_t)hi << 32) | lo;
#else
    register uint64_t imiss;

    asm volatile (
            "csrr %[imiss], mhpmcounter3   \n"
            : [imiss] "=r"(imiss)
            :
            : "memory");

    return imiss;
#endif
#else
    /**********************************************************************************
     *                                      UNSUPPORTED
     *********************************************************************************/
    return 0;

#endif
}

uint64_t cpu_dmiss()
{
#if defined(__x86_64__) || defined(_M_X64)
	return 0;

#elif defined(__riscv)
#if (__SIZEOF_LONG__ == 4)
    register uint32_t lo;
    register uint32_t hi;

    asm volatile (
            "csrr %[lo], mhpmcounter4       \n"
            "csrr %[hi], mhpmcounter4h      \n"
            : [lo] "=r"(lo), [hi] "=r"(hi)
            :
            : "memory");

    /* as there is one cycle between the reading of the lower bits and
     * higher bits, decrement hi if lo was about to wrap */
    if (lo == 0xffffffffUL)
        hi--;

    return ((uint64_t)hi << 32) | lo;
#else
    register uint64_t dmiss;

    asm volatile (
            "csrr %[dmiss], mhpmcounter4   \n"
            : [dmiss] "=r"(dmiss)
            :
            : "memory");

    return dmiss;
#endif
#else
    /**********************************************************************************
     *                                      UNSUPPORTED
     *********************************************************************************/
    return 0;

#endif
}

uint64_t cpu_instructions()
{
#if defined(__x86_64__) || defined(_M_X64)
	return 0;

#elif defined(__riscv)

#if (__SIZEOF_LONG__ == 4)
    register uint32_t lo;
    register uint32_t hi;

    asm volatile (
            "csrr %[lo], minstret                \n"
            "csrr %[hi], minstreth               \n"
            : [lo] "=r"(lo), [hi] "=r"(hi)
            :
            : "memory");

    /* as there is one cycle between the reading of the lower bits and
     * higher bits, decrement hi if lo was about to wrap */
    if (lo == 0xffffffffUL)
        hi--;

    return ((uint64_t)hi << 32) | lo;
#else
    register uint64_t instrs;

    asm volatile (
            "csrr %[instrs], minstret          \n"
            : [instrs] "=r"(instrs)
            :
            : "memory");

    return instrs;
#endif
#else
    /**********************************************************************************
     *                                      UNSUPPORTED
     *********************************************************************************/
    return 0;

#endif
}
