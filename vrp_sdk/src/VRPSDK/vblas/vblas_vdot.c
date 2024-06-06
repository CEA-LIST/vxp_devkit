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
 *  @file        vblas_vdot.c
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
 *  @func   __vdot_1
 *  @brief  Vector-Vector scalar multiplication (dot product)
 */
static inline void __vdot_1(uintptr_t          *x,
                            uintptr_t          *y,
                            const unsigned int ec,
                            const unsigned int evx,
                            const unsigned int evy,
                            int *n,
                            uintptr_t x_bytes,
                            uintptr_t y_bytes)
{
    int _n = *n;
    uintptr_t _x = *x;
    uintptr_t _y = *y;

    while (_n >= 1) {
        ple(P0, _x, 0, evx);
        ple(P1, _y, 0, evy);
        pmul(P0, P0, P1, ec);

        _x += 1*x_bytes;
        _y += 1*y_bytes;
        _n -= 1;

        padd(P31, P0, P31, ec);
    }

    *x     = _x;
    *y     = _y;
    *n     = _n;
}

/**
 *  @func   __vdot_4
 *  @brief  Vector-Vector scalar multiplication (dot product)
 */
static inline void __vdot_4(uintptr_t          *x,
                            uintptr_t          *y,
                            const unsigned int ec,
                            const unsigned int evx,
                            const unsigned int evy,
                            int *n,
                            uintptr_t x_bytes,
                            uintptr_t y_bytes)
{
    int _n = *n;
    uintptr_t _x = *x;
    uintptr_t _y = *y;

    while (_n >= 4) {
        ple(P0, _x, 0, evx);
        ple(P4, _y, 0, evy);
        ple(P1, _x, 1, evx);
        ple(P5, _y, 1, evy);
        ple(P2, _x, 2, evx);
        ple(P6, _y, 2, evy);
        ple(P3, _x, 3, evx);
        ple(P7, _y, 3, evy);
        pmul(P0, P0, P4, ec);
        pmul(P1, P1, P5, ec);
        pmul(P2, P2, P6, ec);
        padd(P0, P0, P31, ec);
        pmul(P3, P3, P7, ec);
        padd(P1, P0, P1, ec);
        padd(P2, P1, P2, ec);

        _x += 4*x_bytes;
        _y += 4*y_bytes;
        _n -= 4;

        padd(P31, P2, P3, ec);
    }

    *x = _x;
    *y = _y;
    *n = _n;
}


/**
 *  @func   __vdot_8_iter
 *  @brief  Vector-Vector scalar multiplication (dot product)
 */
static inline __ALWAYS_INLINE__ void __vdot_8_iter(
                             uintptr_t          *x,
                             uintptr_t          *y,
                             const unsigned int ec,
                             const unsigned int evx,
                             const unsigned int evy,
                             int *n,
                             const uintptr_t x_bytes,
                             const uintptr_t y_bytes,
                             int preamble, int postamble)
{
    int _n = *n;
    uintptr_t _x = *x;
    uintptr_t _y = *y;
    uintptr_t x_incr = 8*x_bytes;
    uintptr_t y_incr = 8*y_bytes;
    if(vblas_likely(_n >= 8)){
        if(preamble) {
            ple(P0, _x, 0, evx);
            ple(P8, _y, 0, evy);
            pcvt_d_p(P24, 0.0);
            ple(P1, _x, 1, evx);
            ple(P9, _y, 1, evy);
            pcvt_d_p(P25, 0.0);
            ple(P2, _x, 2, evx);
            ple(P10, _y, 2, evy);
            pcvt_d_p(P26, 0.0);
            ple(P3, _x, 3, evx);
            ple(P11, _y, 3, evy);
            pcvt_d_p(P27, 0.0);
            ple(P4, _x, 4, evx);
            ple(P12, _y, 4, evy);
            pcvt_d_p(P28, 0.0);
            ple(P5, _x, 5, evx);
            ple(P13, _y, 5, evy);
            pcvt_d_p(P29, 0.0);
            ple(P6, _x, 6, evx);
            ple(P14, _y, 6, evy);
            pcvt_d_p(P30, 0.0);
            // P31 already initialized

            ple(P7, _x, 7, evx);
            _x += x_incr;
            ple(P15, _y, 7, evy);
            _y += y_incr;            
        }
        if(vblas_likely(_n >= 2*8)) {
            if(preamble){
                pmul(P16, P0, P8, ec);
                ple(P0, _x, 0, evx);
                ple(P8, _y, 0, evy);
                pmul(P17, P1, P9, ec);
                ple(P1, _x, 1, evx);
                ple(P9, _y, 1, evy);
                pmul(P18, P2, P10, ec);
                ple(P2, _x, 2, evx);
                ple(P10, _y, 2, evy);
                pmul(P19, P3, P11, ec);
                ple(P3, _x, 3, evx);
                ple(P11, _y, 3, evy);
                pmul(P20, P4, P12, ec);
                ple(P4, _x, 4, evx);
                ple(P12, _y, 4, evy);
                pmul(P21, P5, P13, ec);
                ple(P5, _x, 5, evx);
                ple(P13, _y, 5, evy);
                pmul(P22, P6, P14, ec);
                ple(P6, _x, 6, evx);
                ple(P14, _y, 6, evy);
                pmul(P23, P7, P15, ec);
                ple(P7, _x, 7, evx);
                _x += x_incr;
                ple(P15, _y, 7, evy);
                _y += y_incr;    
            }
            while (vblas_likely((!postamble && (_n >= 8)) || (_n >= 3*8))) {
                padd(P24, P24, P16, ec);
                pmul(P16, P0, P8, ec);
                ple(P0, _x, 0, evx);
                ple(P8, _y, 0, evy);
                padd(P25, P25, P17, ec);
                pmul(P17, P1, P9, ec);
                ple(P1, _x, 1, evx);
                ple(P9, _y, 1, evy);
                padd(P26, P26, P18, ec);
                pmul(P18, P2, P10, ec);
                ple(P2, _x, 2, evx);
                ple(P10, _y, 2, evy);
                padd(P27, P27, P19, ec);
                pmul(P19, P3, P11, ec);
                ple(P3, _x, 3, evx);
                ple(P11, _y, 3, evy);
                padd(P28, P28, P20, ec);
                pmul(P20, P4, P12, ec);
                ple(P4, _x, 4, evx);
                ple(P12, _y, 4, evy);
                padd(P29, P29, P21, ec);
                pmul(P21, P5, P13, ec);
                ple(P5, _x, 5, evx);
                ple(P13, _y, 5, evy);
                padd(P30, P30, P22, ec);
                pmul(P22, P6, P14, ec);
                ple(P6, _x, 6, evx);
                ple(P14, _y, 6, evy);
                padd(P31, P31, P23, ec);
                _n -= 8;
                pmul(P23, P7, P15, ec);
                ple(P7, _x, 7, evx);
                _x += x_incr;
                ple(P15, _y, 7, evy);
                _y += y_incr;   
            }
            if(postamble){
                padd(P24, P24, P16, ec);
                padd(P25, P25, P17, ec);
                padd(P26, P26, P18, ec);
                padd(P27, P27, P19, ec);
                padd(P28, P28, P20, ec);
                padd(P29, P29, P21, ec);
                padd(P30, P30, P22, ec);
                _n -= 8;               
                padd(P31, P31, P23, ec);
            }
        }
        if(postamble){
            pmul(P16, P0, P8, ec);
            pmul(P17, P1, P9, ec);
            pmul(P18, P2, P10, ec);
            pmul(P19, P3, P11, ec);
            pmul(P20, P4, P12, ec);
            pmul(P21, P5, P13, ec);
            pmul(P22, P6, P14, ec);
            pmul(P23, P7, P15, ec);

            padd(P24, P24, P16, ec);
            padd(P25, P25, P17, ec);
            padd(P26, P26, P18, ec);
            padd(P27, P27, P19, ec);
            padd(P28, P28, P20, ec);
            padd(P29, P29, P21, ec);
            padd(P30, P30, P22, ec);
            padd(P31, P31, P23, ec);

            _n -= 8; 

            padd(P25, P24, P25, ec);
            padd(P27, P26, P27, ec);
            padd(P29, P28, P29, ec);
            padd(P31, P30, P31, ec);
            padd(P27, P25, P27, ec);
            padd(P31, P29, P31, ec);
            padd(P31, P27, P31, ec);
        }
    }

    *x = _x;
    *y = _y;
    *n = _n;
}

void vdot(int precision, int n,
           const void *x, vpfloat_evp_t x_evp,
           const void *y, vpfloat_evp_t y_evp,
           void *r, vpfloat_evp_t r_evp, char enable_prefetch)
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
    pser_evp(pack_evp(r_evp.bis, vpfloat_get_rm_mem(), r_evp.es, r_evp.stride),
             EVP0);
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(), x_evp.es, x_evp.stride),
             EVP1);
    pser_evp(pack_evp(y_evp.bis, vpfloat_get_rm_mem(), y_evp.es, y_evp.stride),
             EVP2);

    uintptr_t _x = (uintptr_t) x;
    uintptr_t _y = (uintptr_t) y;
    uintptr_t _x_prefetch = (uintptr_t) x;
    uintptr_t _y_prefetch = (uintptr_t) y;
    const uintptr_t x_bytes = VPFLOAT_SIZEOF(x_evp);
    const uintptr_t y_bytes = VPFLOAT_SIZEOF(y_evp);

    const int PREFETCH_ELEM_PER_BLOCK = 128; // must be multiple of 8
    const int PREFETCH_BLOCKS = n/PREFETCH_ELEM_PER_BLOCK;
    const int PREFETCH_REMAINING = n - (PREFETCH_ELEM_PER_BLOCK*PREFETCH_BLOCKS);

#if VBLAS_ENABLE_HWPF /* The above code is enabling hardware performance counters for VBLAS operations. */

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

        // We can load at most 512b over 4 cycles (the core loop is 4 instr)
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

        //printf("wait x: %d / wait y: %d \n", hwpf_throttle_x.nwait, hwpf_throttle_y.nwait);

        int enable_hwpf = (n > 32) ; // if enough data to prefetch
        // TODO : properly handle non aligned vectors (if not aligned to cache line)

        if (enable_hwpf) {
    #if 0
            printf("Enabling HW prefetch: stride=%d / nblocks=%d / nlines=%d\n",
                    hwpf_param_x.stride,
                    hwpf_param_x.nblocks,
                    hwpf_param_x.nlines);
    #endif

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

            const uintptr_t _addr_x = _x_prefetch;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
                
            const uintptr_t _addr_y = _y_prefetch;
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



    //  Initialize to 0 the result accumulator
    const double zero = 0.0;
    pcvt_d_p(P31, zero);

    for(int i=0; i<PREFETCH_BLOCKS-1; i++, n-=PREFETCH_ELEM_PER_BLOCK) {
#if VBLAS_ENABLE_HWPF
        if ( enable_prefetch ) {
            //  Prefetch next block
            const uintptr_t _addr_x = _x_prefetch + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);
            _x_prefetch = _addr_x;

            const uintptr_t _addr_y = _y_prefetch + PREFETCH_ELEM_PER_BLOCK*y_bytes;
            hwpf_wait(HWPF_ENGINE_Y_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_Y_ID, _addr_y, HWPF_ENABLE);
            _y_prefetch = _addr_y;
        }
#endif        
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(i==0){
            // First block, preamble to initialize the pipeline,
            //   no postamble as we know there is one more block
            __vdot_8_iter(&_x, &_y, EC0, EVP1, EVP2,
                          &elems, x_bytes, y_bytes,
                          1, 0);
        } else {
            // steady state, no preamble nor postamble
            __vdot_8_iter(&_x, &_y, EC0, EVP1, EVP2,
                          &elems, x_bytes, y_bytes,
                          0, 0);
        }            
    }
    if(PREFETCH_BLOCKS > 0) {
#if VBLAS_ENABLE_HWPF
        if ( enable_prefetch ) {
            //  Prefetch last block
            const uintptr_t _addr_x = _x_prefetch + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
                
            const uintptr_t _addr_y = _y_prefetch + PREFETCH_ELEM_PER_BLOCK*y_bytes;
            hwpf_wait(HWPF_ENGINE_Y_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_Y_REM_ID, _addr_y, HWPF_ENABLE);
        }
#endif  
        int elems = PREFETCH_ELEM_PER_BLOCK;
        if(PREFETCH_BLOCKS > 1){
            // Work has been done before, no preamble, 
            //  last iteration on prefetch blocks
            __vdot_8_iter(&_x, &_y, EC0, EVP1, EVP2,
                          &elems, x_bytes, y_bytes,
                          0, 1);
        } else {
            // There is a single block, we need to do the full procedure
            __vdot_8_iter(&_x, &_y, EC0, EVP1, EVP2,
                          &elems, x_bytes, y_bytes,
                          1, 1);
        }    
        n-=PREFETCH_ELEM_PER_BLOCK;
    }
    __vdot_8_iter(&_x, &_y, EC0, EVP1, EVP2,
                  &n, x_bytes, y_bytes,
                  1, 1);
    if (n >= 4) {
        __vdot_4(&_x, &_y, EC0, EVP1, EVP2,
                 &n, x_bytes, y_bytes);
    }
    if (n >= 1) {
        __vdot_1(&_x, &_y, EC0, EVP1, EVP2,
                 &n, x_bytes, y_bytes);
    }

    /* Forward back the result through memory */
    pse(P31, (uintptr_t)r, 0, EVP0);

    /* restore environment */
    pser_ec(old_ec0, EC0);
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);
    pser_evp(old_evp2, EVP2);

    VBLASPERFMONITOR_FUNCTION_END;
}

float sdot(int n, const float *x, vpfloat_off_t x_inc,
           const float *y, vpfloat_off_t y_inc, char enable_prefetch)
{
    float ret;

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_FLOAT;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_FLOAT;
    y_env.stride = y_inc;

    vdot(32, n,
         (void*)x, x_env,
         (void*)y, y_env,
         (void*)&ret, VPFLOAT_EVP_FLOAT, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;

    return ret;
}

double ddot(int n, const double *x, vpfloat_off_t x_inc,
            const double *y, vpfloat_off_t y_inc, char enable_prefetch)
{
    double ret;

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vpfloat_evp_t x_env = VPFLOAT_EVP_DOUBLE;
    x_env.stride = x_inc;

    vpfloat_evp_t y_env = VPFLOAT_EVP_DOUBLE;
    y_env.stride = y_inc;

    vdot(64, n,
         (void*)x, x_env,
         (void*)y, y_env,
         (void*)&ret, VPFLOAT_EVP_DOUBLE, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;

    return ret;
}
