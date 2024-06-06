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
 *  @file        vmath.c
 *  @author      Cesar Fuguet-Tortolero
 */
#include <stdlib.h>
#include <string.h>
// #if DEBUG
#include <stdio.h>
// #endif
#include "VRPSDK/vmath.h"
#include "VRPSDK/vutils.h"

vpfloat_rm_t __vpfloat_rm_mem  = RNE; /* memory rounding mode */
vpfloat_rm_t __vpfloat_rm_comp = RNE; /* compute rounding mode*/

uint64_t bias(vpfloat_es_t exponent_size);


void* vadd(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env)
{
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_ec0 = pger_ec(EC0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));

    //  r = x + y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    padd(P0, P0, P1, EC0);

    //  write-back the result
    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);
    pse(P0, (uintptr_t)r, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_ec(old_ec0, EC0);

    return r;
}

void* vsub(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env)
{
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_ec0 = pger_ec(EC0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));

    //  r = x - y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    psub(P0, P0, P1, EC0);

    //  write-back the result
    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);
    pse(P0, (uintptr_t)r, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_ec(old_ec0, EC0);

    return r;
}

void* vmul(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env)
{
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_ec0 = pger_ec(EC0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));

    //  r = x * y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    pmul(P0, P0, P1, EC0);

    //  write-back the result
    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);
    pse(P0, (uintptr_t)r, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_ec(old_ec0, EC0);

    return r;
}

void* vreciprocal(int precision,
                  const void *x, vpfloat_evp_t x_env,
                  void *r, vpfloat_evp_t r_env)
{
    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_ec0 = pger_ec(EC0);

    //  compute the initial approximation
    const double x0 = 1.0 / vtod(x, x_env);

    //  compute 1/x using the Newton-Raphson method
    //    x(k+1) = x(k) * (2 - (d * x(k)))
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);

    ple(P0, (uintptr_t)x, 0, EVP0);
    pcvt_d_p(P1, dtoraw(x0));
    pcvt_d_p(P2, dtoraw(2.0));
    for (int i = 0; i < 9; i++) {
        pmul(P3, P0, P1, EC0); // d * x(k)
        psub(P3, P2, P3, EC0); // 2 - (d * x(k))
        pmul(P1, P1, P3, EC0); // x(k) * (2 - (d * x(k)))
    }

    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);
    pse(P1, (uintptr_t)r, 0, EVP0);

#if DEBUG
    {
        double _x = vtod(x, x_env);
        double _r = vtod(r, r_env);
        printf("\t(debug) 1 / x = 1 / %e = %e\n", _x, _r);
    }
#endif

    //  restore environment
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);

    return r;
}

void* vdiv(int precision,
           const void *x, vpfloat_evp_t x_env,
           const void *y, vpfloat_evp_t y_env,
           void *r, vpfloat_evp_t r_env)
{
    //  r = x * (1 / y)
    void *t = vreciprocal(precision, y, y_env, NULL, r_env); 
    r = vmul(precision, x, x_env, t, r_env, r, r_env);

#if DEBUG
    {
        printf("\t(debug) x / y = %e / %e = %e\n",
               vtod(x, x_env),
               vtod(y, y_env),
               vtod(r, r_env));
    }
#endif

    free(t);

    return r;
}

unsigned int vsqrt_segments(int p) {
  unsigned int iters=0;
  if (p<6)       iters=6;
  else if (p<20)  iters=7;
  else if (p<40)  iters=8;
  else if (p<100) iters=9;
  else if (p<180) iters=10;
  else if (p<360) iters=11;
  else if (p<400) iters=13;
  else iters=14;
  return(iters);
}

void* vsqrt(int precision,
           const void *x, vpfloat_evp_t x_env,
           void *r, vpfloat_evp_t r_env)
{
    const uint64_t old_evp0 = pger_evp(EVP0);  
    const uint64_t old_evp1 = pger_evp(EVP1);  
    const uint64_t old_ec0 = pger_ec(EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) {
        r = (void*)malloc(VPFLOAT_SIZEOF(r_env));
        memset(r, 0, VPFLOAT_SIZEOF(r_env));
    };

    // initialize environments
    pser_evp(pack_evp(x_env.bis, vpfloat_get_rm_comp(), x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(r_env.bis, vpfloat_get_rm_comp(), r_env.es, r_env.stride), EVP1);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);
    
    ple(P0, (uintptr_t)x, 0, EVP0); // P0 = x;
  
    short is_odd; int c;
    int expo_a = pmv_p_x_exp(P0);

    if ((expo_a & 1)==0) { // even
        is_odd=0;
    } else {
        is_odd=1;
    }

    c=(expo_a + is_odd)/2; // ceil(exp_a/2)

    pmv_x_p_exp(P0, expo_a-2*c);
    
    pcvt_d_p(P1, dtoraw(1.0)); // y = 1.0
    pcvt_d_p(P2, dtoraw(3.0)); // constante 3.0

    for (unsigned int m=0; m<vsqrt_segments(precision); m++) {
        // printf(".... loop start %d -> %s - %lf\n", m, vregtostring(P1, EVP0), pcvt_p_d(P1, EVP0));

        // y=y*(VPFloat(3.0)-(xprime*y*y));
        pmul(P4, P1, P1, EC0); // tmp = y * y

        pmul(P5, P0, P4, EC0); // tmp = xprime * y * y

        psub(P4, P2,  P5, EC0); // tmp = 3.0 - tmp

        pmul(P1, P1, P4, EC0); // y = y * tmp
  
        // div y by 2
        pmv_x_p_exp(P1,  pmv_p_x_exp(P1) - 1);

    }

    // y.exponent(y.exponent()+c);
    pmv_x_p_exp(P1, pmv_p_x_exp(P1) + c);

    // y=xprime*y;
    pmul(P4, P0, P1, EC0);


    pse(P4, (uintptr_t)r, 0, EVP1);

    //  restore environment
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp1, EVP1);
    pser_evp(old_evp0, EVP0);
    
    return r;
}

void* vabs(int precision, const void *x, vpfloat_evp_t x_env,
             void *r, vpfloat_evp_t r_env)
{
    const double CONST_MINUS_1 = -1.0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);  
    const uint64_t old_evp2 = pger_evp(EVP2);  
    const uint64_t old_ec0 = pger_ec(EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) {
        r = (void*)malloc(VPFLOAT_SIZEOF(r_env));
        memset(r, 0, VPFLOAT_SIZEOF(r_env));
    };

    //  initialize environments
    pser_evp(pack_evp(x_env.bis, vpfloat_get_rm_comp(), x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(r_env.bis, vpfloat_get_rm_comp(), r_env.es, r_env.stride), EVP1);
    pser_evp(pack_evp(VPFLOAT_EVP_DOUBLE.bis, vpfloat_get_rm_comp(), VPFLOAT_EVP_DOUBLE.es, VPFLOAT_EVP_DOUBLE.stride), EVP2);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);


    // load x into p register
    ple(P0, (uintptr_t)x, 0, EVP0);

    // Retrieve the sign of the vpfloat scalar
    const uint64_t l_sign = pmv_p_x_sign(P0);

    if ( l_sign == 0 ) {
        // it is positive so y=|x|=x
        pse(P0, (uintptr_t)r, 0, EVP1);
    } else {
        // it is negative so y=|x|=-x

        // Load -1 in a p register
        
        ple(P1, (uintptr_t)&CONST_MINUS_1, 0, EVP2);

        // x = x * -1 = -x
        pmul(P0, P0, P1, EC0);

        // y = -x
        pse(P0, (uintptr_t)r, 0, EVP1);
    }

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_evp(old_evp2, EVP2);
    pser_ec(old_ec0, EC0);

    return r;
}

void vassign(const void *x, vpfloat_evp_t x_env,
             void *y, vpfloat_evp_t y_env)
{
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  y = x
    ple(P0, (uintptr_t)x, 0, EVP0);
    pse(P0, (uintptr_t)y, 0, EVP1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
}

void * vinvertsign(int precision, const void *x, vpfloat_evp_t x_env)
{
    double l_minus_1 = -1.0;

    return vmul(precision,
                x, x_env, 
                (void *)&l_minus_1, VPFLOAT_EVP_DOUBLE,
                NULL, x_env);
}

/*
 * Return the integer representing one chunk in the
 * mantissa of one variable-precision scalar
 */
uint64_t vmantissaChunk(const uint16_t chunk, const void *x, vpfloat_evp_t x_env) {
    uint64_t l_mantissa_chunk = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);

    ple(P0, (uintptr_t)x, 0, EVP0);

    switch(chunk) {
        case P_CHUNK0:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK0, P0);
            break;
        case P_CHUNK1:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK1, P0);
            break;
        case P_CHUNK2:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK2, P0);
            break;
        case P_CHUNK3:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK3, P0);
            break;
        case P_CHUNK4:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK4, P0);
            break;
        case P_CHUNK5:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK5, P0);
            break;
        case P_CHUNK6:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK6, P0);
            break;            
        case P_CHUNK7:
            l_mantissa_chunk = pmv_p_x_chunk(P_CHUNK7, P0);
            break;            
        default:
            printf("%s: chunk %d is not valid.\n", __func__, chunk);
            exit(1);
            break;            
    }

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return l_mantissa_chunk;
}

/*
 * Set the integer value of one chunk in the mantissa of one variable-precision scalar
 */
void vmantissaChunk_set(const uint16_t chunk, const void *x, vpfloat_evp_t x_env, const uint64_t value) {
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);

    ple(P0, (uintptr_t)x, 0, EVP0);

    switch(chunk) {
        case P_CHUNK0:
            pmv_x_p_chunk(P0, P_CHUNK0, value);
            break;
        case P_CHUNK1:
            pmv_x_p_chunk(P0, P_CHUNK1, value);
            break;
        case P_CHUNK2:
            pmv_x_p_chunk(P0, P_CHUNK2, value);
            break;
        case P_CHUNK3:
            pmv_x_p_chunk(P0, P_CHUNK3, value);
            break;
        case P_CHUNK4:
            pmv_x_p_chunk(P0, P_CHUNK4, value);
            break;
        case P_CHUNK5:
            pmv_x_p_chunk(P0, P_CHUNK5, value);
            break;
        case P_CHUNK6:
            pmv_x_p_chunk(P0, P_CHUNK6, value);
            break;
        case P_CHUNK7:
            pmv_x_p_chunk(P0, P_CHUNK7, value);
            break;
        default:
            printf("%s: chunk %d is not valid.\n", __func__, chunk);
            exit(1);
            break;            
    }

    pse(P0, (uintptr_t)x, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
}

int64_t vexponent(const void *x, vpfloat_evp_t x_env) {
    uint64_t l_exponent = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);

    ple(P0, (uintptr_t)x, 0, EVP0);

    l_exponent = pmv_p_x_exp(P0);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return (int64_t)l_exponent;
}


void vexponent_set(const void *x, vpfloat_evp_t x_env, const uint64_t exponent) {
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);

    ple(P0, (uintptr_t)x, 0, EVP0);
    pmv_x_p_exp(P0, exponent);
    pse(P0, (uintptr_t)x, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
}


void pow2(int n, void * r, vpfloat_evp_t r_env) {
    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    // Initialize chunk 0 only an set number of chunk (len) to 1
    pmv_x_p_chunk(P0, P_CHUNK0, 0);
    pmv_x_p_len(P0, 0); // Set len to 1 chunk by setting it to 0    
    pmv_x_p_exp(P0, n); // Set exponent to n
    pmv_x_p_sign(P0, 0);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);

    //  write-back the result
    pse(P0, (uintptr_t)r, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);
}


uint64_t vcompare_gt(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_gt(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

uint64_t vcompare_geq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_geq(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

uint64_t vcompare_lt(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_lt(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

uint64_t vcompare_leq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_leq(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

uint64_t vcompare_eq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_eq(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

uint64_t vcompare_neq(const void *x, vpfloat_evp_t x_env,
             const void *y, vpfloat_evp_t y_env) {
    uint64_t l_comparison_result = 0;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(x_env.bis, RM, x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(y_env.bis, RM, y_env.es, y_env.stride), EVP1);

    //  comparison_result = x > y
    ple(P0, (uintptr_t)x, 0, EVP0);
    ple(P1, (uintptr_t)y, 0, EVP1);
    l_comparison_result = pcmp_neq(P0, P1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    return l_comparison_result;
}

void vnorm( int precision, const void *x, vpfloat_evp_t x_env,
            void *r, vpfloat_evp_t r_env) {

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);  
    const uint64_t old_ec0 = pger_ec(EC0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) {
        r = (void*)malloc(VPFLOAT_SIZEOF(r_env));
        memset(r, 0, VPFLOAT_SIZEOF(r_env));
    };

    //  initialize environments
    pser_evp(pack_evp(x_env.bis, vpfloat_get_rm_comp(), x_env.es, x_env.stride), EVP0);
    pser_evp(pack_evp(r_env.bis, vpfloat_get_rm_comp(), r_env.es, r_env.stride), EVP1);
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()), EC0);


    // load x into p register
    ple(P0, (uintptr_t)x, 0, EVP0); // P0 = x.real
    pmul(P1, P0, P0, EVP0); // P1 = x.real * x.real
    ple(P0, (uintptr_t)(x + VPFLOAT_SIZEOF(x_env)), 0, EVP0); // P0 = x.imag
    pmul(P0, P0, P0, EVP0); // P0 = x.imag * x.imag
    padd(P0, P0, P1, EC0); // P0 = x.imag * x.imag + x.real * x.real

    pse(P0, (uintptr_t)r, 0, EVP1);

    //  restore environment
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_ec(old_ec0, EC0);

    vsqrt(precision, r, r_env, r, r_env);
}