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
 *  @file        vblas_vaxpy.c
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
 *  @func   __vaxpy_1
 *  @brief  Vector scaling by alpha scalar and vector addition for one
 *          element
 *  @note   P0 contains the alpha scaling factor
 */
static inline void __vaxpy_1(const uintptr_t    x,
                             const uintptr_t    y,
                             const unsigned int ec,
                             const unsigned int evx,
                             const unsigned int evy)
{
    ple(P1, x, 0, evx);
    ple(P2, y, 0, evy);
    pmul(P1, P1, P0, ec);
    padd(P1, P1, P2, ec);
    pse(P1, y, 0, evy);
}

/**
 *  @func   __vaxpy_4
 *  @brief  Vector scaling by alpha scalar and vector addition for four
 *          contiguous elements
 *  @note   P0 contains the alpha scaling factor
 */
static inline void __vaxpy_4(const uintptr_t    x,
                             const uintptr_t    y,
                             const unsigned int ec,
                             const unsigned int evx,
                             const unsigned int evy)
{
    ple(P1, x, 0, evx);
    ple(P2, x, 1, evx);
    ple(P3, x, 2, evx);
    ple(P4, x, 3, evx);
    pmul(P1, P1, P0, ec);
    ple(P5, y, 0, evy);
    ple(P6, y, 1, evy);
    pmul(P2, P2, P0, ec);
    ple(P7, y, 2, evy);
    ple(P8, y, 3, evy);
    pmul(P3, P3, P0, ec);
    padd(P1, P1, P5, ec);
    pmul(P4, P4, P0, ec);
    padd(P2, P2, P6, ec);
    padd(P3, P3, P7, ec);
    pse(P1, y, 0, evy);
    padd(P4, P4, P8, ec);
    pse(P2, y, 1, evy);
    pse(P3, y, 2, evy);
    pse(P4, y, 3, evy);
}

/**
 *  @func   __vaxpy_8_iter
 *  @brief  Vector scaling by alpha scalar and vector addition for more than 16
 *          contiguous elements
 *  @note   P0 contains the alpha scaling factor
 */
static inline __ALWAYS_INLINE__ void __vaxpy_7_iter(uintptr_t    *x,
                                  uintptr_t    *y_store,
                                  uintptr_t    *y_read,
                                  const unsigned int ec,
                                  const unsigned int evx,
                                  const unsigned int evy,
                                  const int x_bytes,
                                  const int y_bytes,
                                  int *n,
                                  int preamble, int postamble)
{
    uintptr_t _x = *x;
    uintptr_t _y_read = *y_read;
    uintptr_t _y_store = *y_store;
    int _n = *n;

    if(vblas_likely(_n >= 7)){
        // We operate on 4 pipeline stages to maximize spacing between 
        // dependent instructions

        int y_ptr_inc = 7*y_bytes;
        int x_ptr_inc = 7*x_bytes;
        if(preamble){
            // Cycle 0
            ple(P1, _x, 0, evx);
            ple(P2, _x, 1, evx);
            ple(P3, _x, 2, evx);
            ple(P4, _x, 3, evx);
            ple(P5, _x, 4, evx);
            ple(P6, _x, 5, evx);
            ple(P7, _x, 6, evx);
            _x += x_ptr_inc;
        }
        if(vblas_likely(_n >= 7*2)){
            if(preamble){
                ple(P16, _y_read, 0, evy); // Cycle 0
                pmul(P8, P1, P0, ec); // Cycle 0
                ple(P1, _x, 0, evx); // Cycle 1
                ple(P17, _y_read, 1, evy); // Cycle 0
                pmul(P9, P2, P0, ec); // Cycle 0
                ple(P2, _x, 1, evx); // Cycle 1
                ple(P18, _y_read, 2, evy); // Cycle 0
                pmul(P10, P3, P0, ec); // Cycle 0
                ple(P3, _x, 2, evx); // Cycle 1
                ple(P19, _y_read, 3, evy); // Cycle 0
                pmul(P11, P4, P0, ec); // Cycle 0
                ple(P4, _x, 3, evx); // Cycle 1
                ple(P20, _y_read, 4, evy); // Cycle 0
                pmul(P12, P5, P0, ec); // Cycle 0
                ple(P5, _x, 4, evx); // Cycle 1
                ple(P21, _y_read, 5, evy); // Cycle 0
                pmul(P13, P6, P0, ec); // Cycle 0
                ple(P6, _x, 5, evx); // Cycle 1
                ple(P22, _y_read, 6, evy); // Cycle 0
                pmul(P14, P7, P0, ec); // Cycle 0
                _y_read += y_ptr_inc;
                ple(P7, _x, 6, evx); // Cycle 1
                _x += x_ptr_inc;
            }

            if(vblas_likely(_n >= 7*3)){
                if(preamble){
                    padd(P24, P8, P16, ec); // Cycle 0
                    ple(P16, _y_read, 0, evy); // Cycle 1
                    pmul(P8, P1, P0, ec); // Cycle 1
                    ple(P1, _x, 0, evx); // Cycle 2
                    padd(P25, P9, P17, ec); // Cycle 0
                    ple(P17, _y_read, 1, evy); // Cycle 1
                    pmul(P9, P2, P0, ec); // Cycle 1
                    ple(P2, _x, 1, evx); // Cycle 2
                    padd(P26, P10, P18, ec); // Cycle 0
                    ple(P18, _y_read, 2, evy); // Cycle 1
                    pmul(P10, P3, P0, ec); // Cycle 1
                    ple(P3, _x, 2, evx); // Cycle 2
                    padd(P27, P11, P19, ec); // Cycle 0
                    ple(P19, _y_read, 3, evy); // Cycle 1
                    pmul(P11, P4, P0, ec); // Cycle 1
                    ple(P4, _x, 3, evx); // Cycle 2
                    padd(P28, P12, P20, ec); // Cycle 0
                    ple(P20, _y_read, 4, evy); // Cycle 1
                    pmul(P12, P5, P0, ec); // Cycle 1
                    ple(P5, _x, 4, evx); // Cycle 2
                    padd(P29, P13, P21, ec); // Cycle 0
                    ple(P21, _y_read, 5, evy); // Cycle 1
                    pmul(P13, P6, P0, ec); // Cycle 1
                    ple(P6, _x, 5, evx); // Cycle 2
                    padd(P30, P14, P22, ec); // Cycle 0
                    ple(P22, _y_read, 6, evy); // Cycle 1
                    _y_read += y_ptr_inc;
                    pmul(P14, P7, P0, ec); // Cycle 1
                    ple(P7, _x, 6, evx); // Cycle 2
                    _x += x_ptr_inc;
                }
                // Iterate as long as there is enough data :
                // - at least 7*4 in the general case
                // - if we don't execute the postamble, we assume we must work 
                //     as much as possible
                while(vblas_likely((!postamble && _n >=7) || (_n >= 7*4))){
                    pse(P24, _y_store, 0, evy); // Cycle k-3
                    padd(P24, P8, P16, ec); // Cycle k-2
                    ple(P16, _y_read, 0, evy); // Cycle k-1
                    pmul(P8, P1, P0, ec); // Cycle k-1
                    ple(P1, _x, 0, evx); // Cycle k
                    pse(P25, _y_store, 1, evy); // Cycle k-3
                    padd(P25, P9, P17, ec); // Cycle k-2
                    ple(P17, _y_read, 1, evy); // Cycle k-1
                    pmul(P9, P2, P0, ec); // Cycle k-1
                    ple(P2, _x, 1, evx); // Cycle k
                    pse(P26, _y_store, 2, evy); // Cycle k-3
                    padd(P26, P10, P18, ec); // Cycle k-2
                    ple(P18, _y_read, 2, evy); // Cycle k-1
                    pmul(P10, P3, P0, ec); // Cycle k-1
                    ple(P3, _x, 2, evx); // Cycle k
                    pse(P27, _y_store, 3, evy); // Cycle k-3
                    padd(P27, P11, P19, ec); // Cycle k-2
                    ple(P19, _y_read, 3, evy); // Cycle k-1
                    pmul(P11, P4, P0, ec); // Cycle k-1
                    ple(P4, _x, 3, evx); // Cycle k
                    pse(P28, _y_store, 4, evy); // Cycle k-3
                    padd(P28, P12, P20, ec); // Cycle k-2
                    ple(P20, _y_read, 4, evy); // Cycle k-1
                    pmul(P12, P5, P0, ec); // Cycle k-1
                    ple(P5, _x, 4, evx); // Cycle k
                    pse(P29, _y_store, 5, evy); // Cycle k-3
                    padd(P29, P13, P21, ec); // Cycle k-2
                    ple(P21, _y_read, 5, evy); // Cycle k-1
                    _n -= 7;
                    pmul(P13, P6, P0, ec); // Cycle k-1
                    ple(P6, _x, 5, evx); // Cycle k
                    pse(P30, _y_store, 6, evy); // Cycle k-3
                    _y_store += y_ptr_inc;
                    padd(P30, P14, P22, ec); // Cycle k-2
                    ple(P22, _y_read, 6, evy); // Cycle k-1
                    _y_read += y_ptr_inc;
                    pmul(P14, P7, P0, ec); // Cycle k-1
                    ple(P7, _x, 6, evx); // Cycle k
                    _x += x_ptr_inc;
                }
                if (postamble) {
                    pse(P24, _y_store, 0, evy); // Cycle N-2
                    pse(P25, _y_store, 1, evy); // Cycle N-2
                    pse(P26, _y_store, 2, evy); // Cycle N-2
                    pse(P27, _y_store, 3, evy); // Cycle N-2
                    pse(P28, _y_store, 4, evy); // Cycle N-2
                    pse(P29, _y_store, 5, evy); // Cycle N-2
                    pse(P30, _y_store, 6, evy); // Cycle N-2
                    _y_store += y_ptr_inc;
                    _n -= 7;
                }
            }
            if (postamble) {
                padd(P24, P8, P16, ec); // Cycle N-1
                padd(P25, P9, P17, ec); // Cycle N-1
                padd(P26, P10, P18, ec); // Cycle N-1
                padd(P27, P11, P19, ec); // Cycle N-1
                padd(P28, P12, P20, ec); // Cycle N-1
                padd(P29, P13, P21, ec); // Cycle N-1
                padd(P30, P14, P22, ec); // Cycle N-1
                pse(P24, _y_store, 0, evy); // Cycle N-1
                pse(P25, _y_store, 1, evy); // Cycle N-1
                pse(P26, _y_store, 2, evy); // Cycle N-1
                pse(P27, _y_store, 3, evy); // Cycle N-1
                pse(P28, _y_store, 4, evy); // Cycle N-1
                pse(P29, _y_store, 5, evy); // Cycle N-1
                pse(P30, _y_store, 6, evy); // Cycle N-1
                _y_store += y_ptr_inc;
                _n -= 7;
            }
        }
        if (postamble) {
            ple(P16, _y_read, 0, evy); // Cycle N
            pmul(P8, P1, P0, ec); // Cycle N
            ple(P17, _y_read, 1, evy); // Cycle N
            pmul(P9, P2, P0, ec); // Cycle N
            ple(P18, _y_read, 2, evy); // Cycle N
            pmul(P10, P3, P0, ec); // Cycle N
            ple(P19, _y_read, 3, evy); // Cycle N
            pmul(P11, P4, P0, ec); // Cycle N
            ple(P20, _y_read, 4, evy); // Cycle N
            pmul(P12, P5, P0, ec); // Cycle N
            ple(P21, _y_read, 5, evy); // Cycle N
            pmul(P13, P6, P0, ec); // Cycle N
            ple(P22, _y_read, 6, evy); // Cycle N
            pmul(P14, P7, P0, ec); // Cycle N
            _y_read += y_ptr_inc;
            padd(P24, P8, P16, ec); // Cycle N
            padd(P25, P9, P17, ec); // Cycle N
            padd(P26, P10, P18, ec); // Cycle N
            padd(P27, P11, P19, ec); // Cycle N
            padd(P28, P12, P20, ec); // Cycle N
            padd(P29, P13, P21, ec); // Cycle N
            padd(P30, P14, P22, ec); // Cycle N
            pse(P24, _y_store, 0, evy); // Cycle N
            pse(P25, _y_store, 1, evy); // Cycle N
            pse(P26, _y_store, 2, evy); // Cycle N
            pse(P27, _y_store, 3, evy); // Cycle N
            pse(P28, _y_store, 4, evy); // Cycle N
            pse(P29, _y_store, 5, evy); // Cycle N
            pse(P30, _y_store, 6, evy); // Cycle N
            _y_store += y_ptr_inc;
            _n -= 7;
        }
    }

    *y_read = _y_read;
    *y_store = _y_store;
    *n = _n;
    *x = _x; 
}

void vaxpy(int precision, int n,
           const void *alpha, vpfloat_evp_t a_evp,
           const void *x, vpfloat_evp_t x_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* save environment */
    const uint64_t old_ec0 = pger_ec(EC0);
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);
    const uint64_t old_evp2 = pger_evp(EVP2);

    /* set compute and memory environments */
    pser_ec(pack_ec(precision, vpfloat_get_rm_comp()),
            EC0);
    pser_evp(pack_evp(a_evp.bis, vpfloat_get_rm_mem(), a_evp.es, a_evp.stride),
             EVP0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(), x_evp.es, x_evp.stride),
             EVP1);
    pser_evp(pack_evp(y_evp.bis, vpfloat_get_rm_mem(), y_evp.es, y_evp.stride),
             EVP2);

    uintptr_t x_ptr = (uintptr_t) x;
    uintptr_t y_read_ptr = (uintptr_t) y;
    uintptr_t y_store_ptr = (uintptr_t) y;
    uintptr_t x_prefetch_ptr = (uintptr_t) x;
    uintptr_t y_prefetch_ptr = (uintptr_t) y;
    int _n = n;
    const uintptr_t x_bytes = VPFLOAT_SIZEOF(x_evp);
    const uintptr_t y_bytes = VPFLOAT_SIZEOF(y_evp);

    ple(P0, (uintptr_t)alpha, 0, EVP0);

    // A block should be a multiple of 64B for good prefetch results
    // As we don't know the size of vp elems, 
    //  we suppose it is at least a multiple of 32b aligned on 32b
    // To maximize performance, we should also have a multiple
    //  of 7 elements as is iterates faster this way
    // So we choose 112 elements per block (divisble by 7 and 
    //  leads to 64B divisible blocks for 32b aligned vpfloats)
    const int PREFETCH_ELEM_PER_BLOCK = 16*7; // WARNING : won't work if not divisible by 7!
    const int PREFETCH_BLOCKS = n/PREFETCH_ELEM_PER_BLOCK;
    const int PREFETCH_REMAINING = n - (PREFETCH_ELEM_PER_BLOCK*PREFETCH_BLOCKS);

#if VBLAS_ENABLE_HWPF
    const int HWPF_ENGINE_X_ID = 0;
    const int HWPF_ENGINE_Y_ID = 1;
    const int HWPF_ENGINE_X_REM_ID = 2;
    const int HWPF_ENGINE_Y_REM_ID = 3;

    hwpf_engine_params_t hwpf_param_x;
    hwpf_engine_params_t hwpf_param_y;
    hwpf_engine_params_t hwpf_param_x_rem;
    hwpf_engine_params_t hwpf_param_y_rem;

    if ( enable_prefetch ) {
        hwpf_param_x.stride = 1;
        hwpf_param_x.nlines = 1;
        hwpf_param_x.nblocks = __bytes_to_cachelines(PREFETCH_ELEM_PER_BLOCK*x_bytes);

        hwpf_param_y.stride = 1;
        hwpf_param_y.nlines = 1;
        hwpf_param_y.nblocks = __bytes_to_cachelines(PREFETCH_ELEM_PER_BLOCK*y_bytes);

        hwpf_param_x_rem.stride = 1;
        hwpf_param_x_rem.nlines = 1;
        hwpf_param_x_rem.nblocks = __bytes_to_cachelines(PREFETCH_REMAINING*x_bytes);

        hwpf_param_y_rem.stride = 1;
        hwpf_param_y_rem.nlines = 1;
        hwpf_param_y_rem.nblocks = __bytes_to_cachelines(PREFETCH_REMAINING*y_bytes);

        // We can load at most (128*5)b over 5 cycles (the core loop is 5 instr)
        // So we should send at most one req every 4 cycles for both prefetchers
        //  ie min total wait=8
        int prefetch_wait = (x_bytes+y_bytes) > 64 ? 4 : (4*64)/(x_bytes+y_bytes);
        // account for our imperfect memory hierarchy
        int prefetch_factor = 4;

        hwpf_engine_throttle_t hwpf_throttle_x;
        hwpf_throttle_x.nwait = prefetch_factor+(prefetch_wait+((prefetch_wait*y_bytes)/x_bytes));
        hwpf_throttle_x.ninflight = 0;
        hwpf_engine_throttle_t hwpf_throttle_y;
        hwpf_throttle_y.nwait = prefetch_factor+(prefetch_wait+((prefetch_wait*x_bytes)/y_bytes));
        hwpf_throttle_y.ninflight = 0;

        int enable_hwpf = (n > 32) ; // if enough data to prefetch
        // TODO : properly handle non aligned vectors (if not aligned to cache line)

        if (enable_hwpf) {
            hwpf_set_params(HWPF_ENGINE_X_ID, &hwpf_param_x);
            hwpf_set_throttle(HWPF_ENGINE_X_ID, &hwpf_throttle_x);
            hwpf_set_params(HWPF_ENGINE_Y_ID, &hwpf_param_y);
            hwpf_set_throttle(HWPF_ENGINE_Y_ID, &hwpf_throttle_y);

            if(PREFETCH_BLOCKS > 0) {
                // We use REM engines to prefetch first block while we use it
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle_x);
                hwpf_set_params(HWPF_ENGINE_Y_REM_ID, &hwpf_param_y);
                hwpf_set_throttle(HWPF_ENGINE_Y_REM_ID, &hwpf_throttle_y);
            } else {
                // Else it is also the last block !
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle_x);
                hwpf_set_params(HWPF_ENGINE_Y_REM_ID, &hwpf_param_y_rem);
                hwpf_set_throttle(HWPF_ENGINE_Y_REM_ID, &hwpf_throttle_y);            
            }

            const uintptr_t _addr_x = x_prefetch_ptr;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
                
            const uintptr_t _addr_y = y_prefetch_ptr;
            hwpf_wait(HWPF_ENGINE_Y_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_Y_REM_ID, _addr_y, HWPF_ENABLE);

            if(PREFETCH_BLOCKS > 0) {
                // If the first block is not also the last one, we prepare prefetch
                // of the last one
                hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
                hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle_x);
                hwpf_set_params(HWPF_ENGINE_Y_REM_ID, &hwpf_param_y_rem);
                hwpf_set_throttle(HWPF_ENGINE_Y_REM_ID, &hwpf_throttle_y);             
            }

        }
    }
#endif

    for(int i=0; i<PREFETCH_BLOCKS-1; i++, _n-=PREFETCH_ELEM_PER_BLOCK) {
#if VBLAS_ENABLE_HWPF
        if ( enable_prefetch ) {
            //  Prefetch next block
            const uintptr_t _addr_x = x_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
            x_prefetch_ptr = _addr_x;
            
            const uintptr_t _addr_y = y_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*y_bytes;
            hwpf_wait(HWPF_ENGINE_Y_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_Y_ID, _addr_y, HWPF_ENABLE);
            y_prefetch_ptr = _addr_y;
        }
#endif        
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(i==0){
            // First block, preamble to initialize the pipeline,
            //   no postamble as we know there is one more block
            __vaxpy_7_iter(&x_ptr, &y_read_ptr, &y_store_ptr,
                   EC0, EVP1, EVP2,
                   x_bytes, y_bytes, &elems, 1, 0);
        } else {
            // steady state, no preamble nor postamble
            __vaxpy_7_iter(&x_ptr, &y_read_ptr, &y_store_ptr,
                   EC0, EVP1, EVP2,
                   x_bytes, y_bytes, &elems, 0, 0);
        }
    }
    if(PREFETCH_BLOCKS > 0) {
#if VBLAS_ENABLE_HWPF
        if ( enable_prefetch ) {
            //  Prefetch last block
            const uintptr_t _addr_x = x_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
            x_prefetch_ptr = _addr_x;
                
            const uintptr_t _addr_y = y_prefetch_ptr + PREFETCH_ELEM_PER_BLOCK*y_bytes;
            hwpf_wait(HWPF_ENGINE_Y_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_Y_REM_ID, _addr_y, HWPF_ENABLE);
            y_prefetch_ptr = _addr_y;
        }
#endif  
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(PREFETCH_BLOCKS > 1){
            // Work has been done before, no preamble, 
            //  last iteration on prefetch blocks
            __vaxpy_7_iter(&x_ptr, &y_read_ptr, &y_store_ptr,
                       EC0, EVP1, EVP2,
                       x_bytes, y_bytes, &elems, 0, 1);
        } else {
            // There is a single block, we need to do the full procedure
            __vaxpy_7_iter(&x_ptr, &y_read_ptr, &y_store_ptr,
                       EC0, EVP1, EVP2,
                       x_bytes, y_bytes, &elems, 1, 1);
        }
        _n-=PREFETCH_ELEM_PER_BLOCK;
    }

    // One more call for last block, we cannot optimize with previous ones
    //  as number of elements is not the same
    __vaxpy_7_iter(&x_ptr, &y_read_ptr, &y_store_ptr,
                   EC0, EVP1, EVP2,
                   x_bytes, y_bytes, &_n, 1 ,1);
    
    if (_n >= 4) {
        __vaxpy_4(x_ptr, y_read_ptr, EC0, EVP1, EVP2);
        _n -= 4;
        x_ptr += 4*x_bytes;
        y_read_ptr += 4*y_bytes;
    }
    while (_n > 0) {
        __vaxpy_1(x_ptr, y_read_ptr, EC0, EVP1, EVP2);
        _n -= 1;
        x_ptr += 1*x_bytes;
        y_read_ptr += 1*y_bytes;
    }
    /* restore environment */
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_evp(old_evp2, EVP2);

    VBLASPERFMONITOR_FUNCTION_END;
}

void saxpy(int n, float alpha,
           const float *x, vpfloat_off_t x_inc,
           float *y, vpfloat_off_t y_inc, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_FLOAT;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_FLOAT;
    y_env.stride = y_inc;

    vaxpy(32, n,
          (void*)&alpha, VPFLOAT_EVP_FLOAT,
          (void*)x, x_env,
          (void*)y, y_env, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;
}

void daxpy(int n, double alpha,
           const double *x, vpfloat_off_t x_inc,
           double *y, vpfloat_off_t y_inc, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_DOUBLE;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_DOUBLE;
    y_env.stride = y_inc;

    vaxpy(64, n,
          (void*)&alpha, VPFLOAT_EVP_DOUBLE,
          (void*)x, x_env,
          (void*)y, y_env, enable_prefetch);

        VBLASPERFMONITOR_FUNCTION_END;
}
