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
 *  @file        vblas_vscal.c
 *  @author      Cesar Fuguet-Tortolero
 */
#include "VRPSDK/vblas.h"
#include "VRPSDK/compiler.h"
#include "VRPSDK/vblas_perfmonitor.h"

inline int __bytes_to_cachelines(int bytes)
{
    return (bytes + BSP_CONFIG_DCACHE_LINE_BYTES - 1)/BSP_CONFIG_DCACHE_LINE_BYTES;
}

/**
 *  @func   __vscal_1
 *  @brief  Vector scaling operation for a single element
 *  @note   P0 contains the alpha scaling factor
 */
static inline void __vscal_1(const uintptr_t    x,
                             const unsigned int ec,
                             const unsigned int evx)
{
    ple(P1, x, 0, evx);
    pmul(P1, P1, P0, ec);
    pse(P1, x, 0, evx);
}

/**
 *  @func   __vscal_4
 *  @brief  Vector scaling operation for four (contiguous) elements
 *  @note   P0 contains the alpha scaling factor
 */
static inline void __vscal_4(const uintptr_t    x,
                             const unsigned int ec,
                             const unsigned int evx)
{
    ple(P1, x, 0, evx);
    ple(P2, x, 1, evx);
    ple(P3, x, 2, evx);
    ple(P4, x, 3, evx);
    pmul(P1, P1, P0, ec);
    pmul(P2, P2, P0, ec);
    pmul(P3, P3, P0, ec);
    pmul(P4, P4, P0, ec);
    pse(P1, x, 0, evx);
    pse(P2, x, 1, evx);
    pse(P3, x, 2, evx);
    pse(P4, x, 3, evx);
}

/**
 *  @func   __vscal_8_iter
 *  @brief  Vector scaling operation for eight (contiguous) elements
 *  @note   P0 contains the alpha scaling factor
 */
static inline __ALWAYS_INLINE__ void __vscal_8_iter(uintptr_t *x_rd_ptr,
                                  uintptr_t *x_wr_ptr,
                                  const unsigned int ec,
                                  const unsigned int evx,
                                  const int x_bytes,
                                  int *n, int preamble, int postamble)
{
    int _n = *n;
    uintptr_t x_rd = *x_rd_ptr;
    uintptr_t x_wr = *x_wr_ptr;
    uintptr_t x_incr = 8*x_bytes;
    // 3-stage pipelining
    if (vblas_likely(_n >= 8)) {
        if(preamble){
            ple(P8,  x_rd, 0, evx);
            ple(P9,  x_rd, 1, evx);
            ple(P10, x_rd, 2, evx);
            ple(P11, x_rd, 3, evx);
            ple(P12, x_rd, 4, evx);
            ple(P13, x_rd, 5, evx);
            ple(P14, x_rd, 6, evx);
            ple(P15, x_rd, 7, evx);
            x_rd  += x_incr;
        }
        if (vblas_likely(_n >= 2*8)){
            if(preamble){
                pmul(P16, P8,  P0, ec);
                ple(P8,  x_rd, 0, evx);
                pmul(P17, P9,  P0, ec);
                ple(P9,  x_rd, 1, evx);
                pmul(P18, P10, P0, ec);
                ple(P10, x_rd, 2, evx);
                pmul(P19, P11, P0, ec);
                ple(P11, x_rd, 3, evx);
                pmul(P20, P12, P0, ec);
                ple(P12, x_rd, 4, evx);
                pmul(P21, P13, P0, ec);
                ple(P13, x_rd, 5, evx);
                pmul(P22, P14, P0, ec);
                ple(P14, x_rd, 6, evx);
                pmul(P23, P15, P0, ec);
                ple(P15, x_rd, 7, evx);
                x_rd  += x_incr;
            }
            while(vblas_likely((!postamble && (_n >= 8)) || (_n >= 3*8))){
                pse(P16,  x_wr,  0, evx);
                pmul(P16, P8,  P0, ec);
                ple(P8,  x_rd, 0, evx);
                pse(P17,  x_wr,  1, evx);
                pmul(P17, P9,  P0, ec);
                ple(P9,  x_rd, 1, evx);
                pse(P18,  x_wr,  2, evx);
                pmul(P18, P10, P0, ec);
                ple(P10, x_rd, 2, evx);
                pse(P19,  x_wr,  3, evx);
                pmul(P19, P11, P0, ec);
                ple(P11, x_rd, 3, evx);
                pse(P20,  x_wr,  4, evx);
                pmul(P20, P12, P0, ec);
                ple(P12, x_rd, 4, evx);
                pse(P21,  x_wr,  5, evx);
                pmul(P21, P13, P0, ec);
                ple(P13, x_rd, 5, evx);
                pse(P22,  x_wr,  6, evx);
                pmul(P22, P14, P0, ec);
                _n -= 8;
                ple(P14, x_rd, 6, evx);
                pse(P23,  x_wr,  7, evx);
                x_wr  += x_incr;
                pmul(P23, P15, P0, ec);
                ple(P15, x_rd, 7, evx);
                x_rd  += x_incr;
            }
            if(postamble){
                pse(P16,  x_wr,  0, evx);
                pse(P17,  x_wr,  1, evx);
                pse(P18,  x_wr,  2, evx);
                pse(P19,  x_wr,  3, evx);
                pse(P20,  x_wr,  4, evx);
                pse(P21,  x_wr,  5, evx);
                pse(P22,  x_wr,  6, evx);
                _n -= 8;
                pse(P23,  x_wr,  7, evx);
                x_wr  += x_incr;
            }
        }
        if(postamble){
                pmul(P16, P8,  P0, ec);
                pmul(P17, P9,  P0, ec);
                pmul(P18, P10, P0, ec);
                pmul(P19, P11, P0, ec);
                pmul(P20, P12, P0, ec);
                pmul(P21, P13, P0, ec);
                pmul(P22, P14, P0, ec);
                _n -= 8;
                pmul(P23, P15, P0, ec);
                pse(P16,  x_wr,  0, evx);
                pse(P17,  x_wr,  1, evx);
                pse(P18,  x_wr,  2, evx);
                pse(P19,  x_wr,  3, evx);
                pse(P20,  x_wr,  4, evx);
                pse(P21,  x_wr,  5, evx);
                pse(P22,  x_wr,  6, evx);
                pse(P23,  x_wr,  7, evx);
                x_wr  += x_incr;
        }
    }
    *n = _n;
    *x_rd_ptr = x_rd;
    *x_wr_ptr = x_wr;
}

void vscal(int precision, int n,
           const void *alpha, vpfloat_evp_t a_evp,
           void *x, vpfloat_evp_t x_evp, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* save environment */
    const uint64_t old_ec0 = pger_ec(EC0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    /* set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()),
            EC0);
    pser_evp(pack_evp(a_evp.bis, vpfloat_get_rm_mem(), a_evp.es, a_evp.stride),
             EVP0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(), x_evp.es, x_evp.stride),
             EVP1);

    uintptr_t x_rd = (uintptr_t)x;
    uintptr_t x_wr = (uintptr_t)x;
    const uintptr_t bytes = VPFLOAT_SIZEOF(x_evp);
    uintptr_t x_prefetch_ptr = (uintptr_t) x;

    // FIXME optimization: return when alpha == 1
    // FIXME optimization: zero the x vector when alpha == 0
    ple(P0, (uintptr_t)alpha, 0, EVP0);

    // A block should be a multiple of 64B for good prefetch results
    // As we don't know the size of vp elems, 
    //  we suppose it is at least a multiple of 32b aligned on 32b
    // To maximize performance, we should also have a multiple
    //  of 8 elements as is iterates faster this way
    // So we choose 128 elements per block (divisble by 8 and 
    //  leads to 64B divisible blocks for any vpfloats)
    const int PREFETCH_ELEM_PER_BLOCK = 16*8; // WARNING : won't work if not divisible by 8!
    const int PREFETCH_BLOCKS = n/PREFETCH_ELEM_PER_BLOCK;
    const int PREFETCH_REMAINING = n - (PREFETCH_ELEM_PER_BLOCK*PREFETCH_BLOCKS);

#if VBLAS_ENABLE_HWPF
    const int HWPF_ENGINE_X_ID = 0;
    const int HWPF_ENGINE_X_REM_ID = 1;

    hwpf_engine_params_t hwpf_param_x;
    hwpf_engine_params_t hwpf_param_x_rem;
    hwpf_engine_throttle_t hwpf_throttle;

    if( enable_prefetch ) {
        hwpf_param_x.stride = 1;
        hwpf_param_x.nlines = 1;
        hwpf_param_x.nblocks = __bytes_to_cachelines(PREFETCH_ELEM_PER_BLOCK*bytes);

        hwpf_param_x_rem.stride = 1;
        hwpf_param_x_rem.nlines = 1;
        hwpf_param_x_rem.nblocks = __bytes_to_cachelines(PREFETCH_REMAINING*bytes);

        hwpf_throttle.nwait = 8;
        hwpf_throttle.ninflight = 0;

        int enable_hwpf = (n > 32) ; // if enough data to prefetch
        // TODO : properly handle non aligned vectors (if not aligned to cache line)

        if (enable_hwpf) {
            hwpf_set_params(HWPF_ENGINE_X_ID, &hwpf_param_x);
            hwpf_set_throttle(HWPF_ENGINE_X_ID, &hwpf_throttle);

            if(PREFETCH_BLOCKS > 0) {
                // We use REM engines to prefetch first block while we use it
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);
            } else {
                // Else it is also the last block !
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);           
            }

            const uintptr_t _addr_x = x_prefetch_ptr;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);

            if(PREFETCH_BLOCKS > 0) {
                // If the first block is not also the last one, we prepare prefetch
                // of the last one
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);    
            }

        }
    }
#endif

    for(int i=0; i<PREFETCH_BLOCKS-1; i++, n-=PREFETCH_ELEM_PER_BLOCK) {
#if VBLAS_ENABLE_HWPF
        if( enable_prefetch ) {
            //  Prefetch next block
            const uintptr_t _addr_x = x_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*bytes;
            hwpf_wait(HWPF_ENGINE_X_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
            x_prefetch_ptr = _addr_x;
        }
#endif        
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(i==0){
            // First block, preamble to initialize the pipeline,
            //   no postamble as we know there is one more block
            __vscal_8_iter(&x_rd,
                           &x_wr,
                           EC0,
                           EVP1,
                           bytes,
                           &elems, 1, 0);
        } else {
            // steady state, no preamble nor postamble
            __vscal_8_iter(&x_rd,
                           &x_wr,
                           EC0,
                           EVP1,
                           bytes,
                           &elems, 0, 0);
        }
    }
    if(PREFETCH_BLOCKS > 0) {
#if VBLAS_ENABLE_HWPF
        if( enable_prefetch ) {
            //  Prefetch last block
            const uintptr_t _addr_x = x_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*bytes;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
            x_prefetch_ptr = _addr_x;
        }
#endif  
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(PREFETCH_BLOCKS > 1){
            // Work has been done before, no preamble, 
            //  last iteration on prefetch blocks
            __vscal_8_iter(&x_rd,
                           &x_wr,
                           EC0,
                           EVP1,
                           bytes,
                           &elems, 0, 1);
        } else {
            // There is a single block, we need to do the full procedure
            __vscal_8_iter(&x_rd,
                           &x_wr,
                           EC0,
                           EVP1,
                           bytes,
                           &elems, 1, 1);
        }
        n-=PREFETCH_ELEM_PER_BLOCK;
    }

    // One more call for last block, we cannot optimize with previous ones
    //  as number of elements is not the same
    __vscal_8_iter(&x_rd,
                   &x_wr,
                   EC0,
                   EVP1,
                   bytes,
                   &n, 1, 1);
    while (n >= 4) {
        __vscal_4(x_rd, EC0, EVP1);
        n    -= 4;
        x_rd += 4*bytes;
    }
    while (n > 0) {
        __vscal_1(x_rd, EC0, EVP1);
        n    -= 1;
        x_rd += 1*bytes;
    }

    /* restore environment */
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    VBLASPERFMONITOR_FUNCTION_END;

}

void sscal(int n, float alpha,
           float *x, vpfloat_off_t x_inc, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_FLOAT;
    x_env.stride = x_inc;

    vscal(32, n,
          (void*)&alpha, VPFLOAT_EVP_FLOAT,
          (void*)x, x_env, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;
}

void dscal(int n, double alpha,
           double *x, vpfloat_off_t x_inc, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_DOUBLE;
    x_env.stride = x_inc;

    vscal(64, n,
          (void*)&alpha, VPFLOAT_EVP_DOUBLE,
          (void*)x, x_env, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;
}
