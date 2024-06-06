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
 *  @file        vutils.h
 *  @author      Cesar Fuguet-Tortolero
 */
#ifndef __VUTILS_H__
#define __VUTILS_H__

#include "vblas.h"
#include "common/bitset.h"
#include "asm/vpfloat.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Convert a double to raw uint 64-bits (no conversion)
 */
static inline uint64_t dtoraw(double x)
{
    union {double d; uint64_t i;} ret;
    ret.d = x;
    return ret.i;
}

/*
 *  Convert a raw 64-bits uint to double (no conversion)
 */
static inline double rawtod(uint64_t x)
{
    union {double d; uint64_t i;} ret;
    ret.i = x;
    return ret.d;
}

/*
 *  Convert a float to raw uint 32-bits (no conversion)
 */
static inline uint32_t ftoraw(float x)
{
    union {float f; uint64_t i;} ret;
    ret.f = x;
    return ret.i;
}

/*
 *  Convert a raw 32-bits uint to float (no conversion)
 */
static inline float rawtof(uint32_t x)
{
    union {float f; uint32_t i;} ret;
    ret.i = x;
    return ret.f;
}

/*
 *  Return the address of element 'index' in the variable-precision vector 'v'
 */
static inline void* vat(const void *v, vpfloat_evp_t env, int index)
{
    intptr_t offset = index*VPFLOAT_SIZEOF(env);
    intptr_t ret    = (intptr_t)v + offset;
    return (void*)ret;
}


/*
 *  Convert variable-precision scalar in memory to double (64 bits)
 */
double vtod(const void *v, vpfloat_evp_t env);

/*
 *  Convert double scalar (64 bits) to variable precision in memory
 */
void dtov(double d, void *v, vpfloat_evp_t env);

/*
 *  Convert variable-precision scalar in memory to float (32 bits)
 */
float vtof(const void *v, vpfloat_evp_t env);

/*
 *  Convert float scalar (32 bits) to variable precision in memory
 */
void ftov(float f, void *v, vpfloat_evp_t env);

/*
 *  Print a variable-precision vector from memory
 */
void vprint_vector(const char* nm, int n,
                   void *v, vpfloat_evp_t env);

/*
 *  Print a variable-precision matrix from memory
 */
void vprint_matrix(const char* nm, int m, int n, int ld,
                   void *mat, vpfloat_evp_t env);

/*
 * Return a string representing the variable-precision scalar stored P register.
 */
char * vregtostring(const unsigned int preg, const unsigned int ev);

/*
 * Return a string representing the variable-precision scalar stored in memory.
 */
char * vtostring(const void *v, vpfloat_evp_t env);

/*
 * Initialize the memory representation of a vpfloat matching
 * environment env and initialize from exponent and 1st mantissa chunk.
 * If memory pointer is NULL, memory for vpfloat data is allocated 
 * by this function.
 */
void * vbuild(const uint64_t mantissa_1st_chunk, const uint64_t exponent, char a_sign, void *r, vpfloat_evp_t r_env);

/**
 * Compute the variable precision transposition of the (n,p) shape matrix a (without copy).
 * Use common/bitset.h
 */
void vtrans(int precision, int n, int p, const void *a, vpfloat_evp_t a_evp);

/**
 * Function returning the number of chunks specified for vpfloat provided as argument.
 */
uint64_t vlen(const void *v, vpfloat_evp_t env);

uint64_t isNaN(const void *v, vpfloat_evp_t env);

uint64_t isInf(const void *v, vpfloat_evp_t env);

#ifdef __cplusplus
}
#endif

#endif /* __VUTILS_H__ */
