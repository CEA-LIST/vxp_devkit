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
 *  @file        vmath.h
 *  @author      Cesar Fuguet-Tortolero
 */
#ifndef __VMATH_H__
#define __VMATH_H__

#include <stdlib.h>
#include "asm/vpfloat.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Definition of rounding modes
 */
extern vpfloat_rm_t __vpfloat_rm_mem; /* memory rounding mode */
extern vpfloat_rm_t __vpfloat_rm_comp; /* compute rounding mode*/

/**
 *  @func   vpfloat_set_rm_mem
 *  @brief  Set the rounding mode for memory operations
 */
volatile static inline void vpfloat_set_rm_mem(vpfloat_rm_t rm)
{
    __vpfloat_rm_mem = rm;
}

/**
 *  @func   vpfloat_set_rm_comp
 *  @brief  Set the rounding mode for compute operations
 */
volatile static inline void vpfloat_set_rm_comp(vpfloat_rm_t rm)
{
    __vpfloat_rm_comp = rm;
}

/**
 *  @func   vpfloat_get_rm_mem
 *  @brief  Get the rounding mode of memory operations
 */
volatile static inline vpfloat_rm_t vpfloat_get_rm_mem()
{
    return __vpfloat_rm_mem;
}

/**
 *  @func   vpfloat_get_rm_comp
 *  @brief  Get the rounding mode of compute operations
 */
volatile static inline vpfloat_rm_t vpfloat_get_rm_comp()
{
    return __vpfloat_rm_comp;
}

/*
 *  Add two variable-precision scalars (r = x + y)
 */
void* vadd(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env);

/*
 *  Substract two variable-precision scalars (r = x - y)
 */
void* vsub(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env);

/*
 *  Multiply two variable-precision scalars (r = x * y)
 */
void* vmul(int precision,
          const void *x, vpfloat_evp_t x_env,
          const void *y, vpfloat_evp_t y_env,
          void *r, vpfloat_evp_t r_env);

/*
 *  Compute the reciprocal of a variable-precision scalar (r = 1 / x)
 */
void* vreciprocal(int precision,
                  const void *x, vpfloat_evp_t x_env,
                  void *r, vpfloat_evp_t r_env);

/*
 *  Divide two variable-precision scalars (r = x / y)
 */
void* vdiv(int precision,
          const void *x, vpfloat_evp_t x_env,
          const void *y, vpfloat_evp_t y_env,
          void *r, vpfloat_evp_t r_env);


/*
 *  Square root of a variable-precision scalar
 */
void* vsqrt(int precision,
            const void *x, vpfloat_evp_t x_env,
            void *r, vpfloat_evp_t r_env);

/*
 *  Copy one variable-precision scalar to another
 */
void vassign(const void *x, vpfloat_evp_t x_env,
             void *y, vpfloat_evp_t y_env);

/*
 * Invert sign of one variable-precision scalar
 */
void * vinvertsign(int precision, const void *x, vpfloat_evp_t x_env);

/*
 * Return the integer representing one chunk in the
 * mantissa of one variable-precision scalar
 */
uint64_t vmantissaChunk(const uint16_t chunk, const void *x, vpfloat_evp_t x_env);

/*
 * Set the integer value of one chunk in the mantissa of one variable-precision scalar
 */
void vmantissaChunk_set(const uint16_t chunk, const void *x, vpfloat_evp_t x_env, const uint64_t value);

/*
 * Return the integer representing the exponent of
 * one variable-precision scalar
 */
int64_t vexponent(const void *x, vpfloat_evp_t x_env);

/*
 * Set the integer representing the exponent of
 * one variable-precision scalar
 */
void vexponent_set(const void *x, vpfloat_evp_t x_env, const uint64_t exponent);

/*
 * Return the value of 2^n in vpfloat format
 */
void pow2(int n, void * r, vpfloat_evp_t x_env);

/*
 * Return 1 if x > y
 */
uint64_t vcompare_gt(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Return 1 if x >= y
 */
uint64_t vcompare_geq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Return 1 if x < y
 */
uint64_t vcompare_lt(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Return 1 if x <= y
 */
uint64_t vcompare_leq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Return 1 if x == y
 */
uint64_t vcompare_eq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Return 1 if x != y
 */
uint64_t vcompare_neq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

/*
 * Build the absolute value of a variable-precision scalar
 */
void* vabs(int precision,
            const void *x, vpfloat_evp_t x_env,
            void *r, vpfloat_evp_t r_env);

uint64_t vcompare_gt(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env);

void vnorm( int precision, const void *x, vpfloat_evp_t x_env,
            void *r, vpfloat_evp_t r_env);

#ifdef __cplusplus
}
#endif

#endif /* __VMATH_H__ */
