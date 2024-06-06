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
 *  @file        vblas_vcopy.c
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
 *  @func   __vcopy_1
 *  @brief  Vector copy operation for four (contiguous) elements
 */
static inline void __vcopy_1(const uintptr_t    x,
                             const uintptr_t    y,
                             const unsigned int evx,
                             const unsigned int evy)
{
    ple(P0, x, 0, evx);
    pse(P0, y, 0, evy);
}

/**
 *  @func   __vcopy_4
 *  @brief  Vector copy operation for four (contiguous) elements
 */
static inline void __vcopy_4(const uintptr_t    x,
                             const uintptr_t    y,
                             const unsigned int evx,
                             const unsigned int evy)
{
    ple(P0, x, 0, evx);
    ple(P1, x, 1, evx);
    ple(P2, x, 2, evx);
    ple(P3, x, 3, evx);
    pse(P0, y, 0, evy);
    pse(P1, y, 1, evy);
    pse(P2, y, 2, evy);
    pse(P3, y, 3, evy);
}

/**
 *  @func   __vcopy_8
 *  @brief  Vector copy operation for eight (contiguous) elements
 */
static inline void __vcopy_8(const uintptr_t    x,
                             const uintptr_t    y,
                             const unsigned int evx,
                             const unsigned int evy)
{
    ple(P0, x, 0, evx);
    ple(P1, x, 1, evx);
    ple(P2, x, 2, evx);
    ple(P3, x, 3, evx);
    ple(P4, x, 4, evx);
    ple(P5, x, 5, evx);
    ple(P6, x, 6, evx);
    ple(P7, x, 7, evx);
    pse(P0, y, 0, evy);
    pse(P1, y, 1, evy);
    pse(P2, y, 2, evy);
    pse(P3, y, 3, evy);
    pse(P4, y, 4, evy);
    pse(P5, y, 5, evy);
    pse(P6, y, 6, evy);
    pse(P7, y, 7, evy);
}

/**
 *  @func   __vcopy_16
 *  @brief  Vector copy operation for sixteen (contiguous) elements
 */
static inline void __vcopy_16(const uintptr_t    x,
                              const uintptr_t    y,
                              const unsigned int evx,
                              const unsigned int evy)
{
    ple(P0, x, 0, evx);
    ple(P1, x, 1, evx);
    ple(P2, x, 2, evx);
    ple(P3, x, 3, evx);
    ple(P4, x, 4, evx);
    ple(P5, x, 5, evx);
    ple(P6, x, 6, evx);
    ple(P7, x, 7, evx);
    pse(P0, y, 0, evy);
    ple(P8, x, 8, evx);
    pse(P1, y, 1, evy);
    ple(P9, x, 9, evx);
    pse(P2, y, 2, evy);
    ple(P10, x, 10, evx);
    pse(P3, y, 3, evy);
    ple(P11, x, 11, evx);
    pse(P4, y, 4, evy);
    ple(P12, x, 12, evx);
    pse(P5, y, 5, evy);
    ple(P13, x, 13, evx);
    pse(P6, y, 6, evy);
    ple(P14, x, 14, evx);
    pse(P7, y, 7, evy);
    ple(P15, x, 15, evx);
    pse(P8, y, 8, evy);
    pse(P9, y, 9, evy);
    pse(P10, y, 10, evy);
    pse(P11, y, 11, evy);
    pse(P12, y, 12, evy);
    pse(P13, y, 13, evy);
    pse(P14, y, 14, evy);
    pse(P15, y, 15, evy);
}

/**
 *  @func   __vcopy_16
 *  @brief  Vector copy operation for sixteen (contiguous) elements
 */
static inline __ALWAYS_INLINE__ void __vcopy_16_iter(uintptr_t *x_ptr,
                                   uintptr_t *y_ptr,
                                   int       *n,
                                   const unsigned int evx,
                                   const unsigned int evy,
                                   const int x_bytes,
                                   const int y_bytes)
{
    uintptr_t x = *x_ptr;
    uintptr_t y = *y_ptr;
    int _n = *n;
    uintptr_t x_incr = 16*x_bytes;
    uintptr_t y_incr = 16*y_bytes;
    if(_n >= 16) {
        ple(P0, x, 0, evx);
        ple(P1, x, 1, evx);
        ple(P2, x, 2, evx);
        ple(P3, x, 3, evx);
        ple(P4, x, 4, evx);
        ple(P5, x, 5, evx);
        ple(P6, x, 6, evx);
        ple(P7, x, 7, evx);
        ple(P8, x, 8, evx);
        ple(P9, x, 9, evx);
        ple(P10, x, 10, evx);
        ple(P11, x, 11, evx);
        ple(P12, x, 12, evx);
        _n -= 16;
        ple(P13, x, 13, evx);
        ple(P14, x, 14, evx);
        ple(P15, x, 15, evx);
        x += x_incr;
        while(_n >= 16) {
            pse(P0, y, 0, evy);
            ple(P0, x, 0, evx);
            pse(P1, y, 1, evy);
            ple(P1, x, 1, evx);
            pse(P2, y, 2, evy);
            ple(P2, x, 2, evx);
            pse(P3, y, 3, evy);
            ple(P3, x, 3, evx);
            pse(P4, y, 4, evy);
            ple(P4, x, 4, evx);
            pse(P5, y, 5, evy);
            ple(P5, x, 5, evx);
            pse(P6, y, 6, evy);
            ple(P6, x, 6, evx);
            pse(P7, y, 7, evy);
            ple(P7, x, 7, evx);
            pse(P8, y, 8, evy);
            ple(P8, x, 8, evx);
            pse(P9, y, 9, evy);
            ple(P9, x, 9, evx);
            pse(P10, y, 10, evy);
            ple(P10, x, 10, evx);
            pse(P11, y, 11, evy);
            ple(P11, x, 11, evx);
            pse(P12, y, 12, evy);
            ple(P12, x, 12, evx);
            _n -= 16;
            pse(P13, y, 13, evy);
            ple(P13, x, 13, evx);
            pse(P14, y, 14, evy);
            ple(P14, x, 14, evx);
            pse(P15, y, 15, evy);
            y += y_incr;
            ple(P15, x, 15, evx);
            x += x_incr;
        }
        pse(P0, y, 0, evy);
        pse(P1, y, 1, evy);
        pse(P2, y, 2, evy);
        pse(P3, y, 3, evy);
        pse(P4, y, 4, evy);
        pse(P5, y, 5, evy);
        pse(P6, y, 6, evy);
        pse(P7, y, 7, evy);
        pse(P8, y, 8, evy);
        pse(P9, y, 9, evy);
        pse(P10, y, 10, evy);
        pse(P11, y, 11, evy);
        pse(P12, y, 12, evy);
        pse(P13, y, 13, evy);
        pse(P14, y, 14, evy);
        pse(P15, y, 15, evy);
        y += y_incr;
    }
    *n = _n;
    *x_ptr = x;
    *y_ptr = y;
}

void vcopy(int n, const void *x, vpfloat_evp_t x_evp,
           void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    /* save environment */
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_evp1 = pger_evp(EVP1);

    /* set memory environments */
    pser_evp(pack_evp(x_evp.bis, vpfloat_get_rm_mem(), x_evp.es, x_evp.stride),
             EVP0);
    pser_evp(pack_evp(y_evp.bis, vpfloat_get_rm_mem(), y_evp.es, y_evp.stride),
             EVP1);

    uintptr_t _x = (uintptr_t)x;
    uintptr_t _y = (uintptr_t)y;
    const uintptr_t x_bytes = VPFLOAT_SIZEOF(x_evp);
    const uintptr_t y_bytes = VPFLOAT_SIZEOF(y_evp);

#if VBLAS_ENABLE_HWPF
    // A block should be a multiple of 64B for good prefetch results
    // As we don't know the size of vp elems,
    //  we suppose it is at least a multiple of 32b aligned on 32b
    // To maximize performance, we should also have a multiple
    //  of 16 elements as is iterates faster this way
    // So we choose 128 elements per block (divisible by 16 and
    //  leads to 64B divisible blocks for any vpfloats)
    const int PREFETCH_ELEM_PER_BLOCK = 16*8; // WARNING : won't work if not divisible by 16!

    if (enable_prefetch && (n > PREFETCH_ELEM_PER_BLOCK)) {
        const int PREFETCH_BLOCKS = n/PREFETCH_ELEM_PER_BLOCK;
        const int PREFETCH_REMAINING = n - (PREFETCH_ELEM_PER_BLOCK*PREFETCH_BLOCKS);

        const int HWPF_ENGINE_X_ID = 0;
        const int HWPF_ENGINE_X_REM_ID = 1;

        hwpf_engine_params_t hwpf_param_x;
        hwpf_engine_params_t hwpf_param_x_rem;
        hwpf_engine_throttle_t hwpf_throttle;

        hwpf_param_x.stride = 1;
        hwpf_param_x.nlines = 1;
        hwpf_param_x.nblocks = __bytes_to_cachelines(PREFETCH_ELEM_PER_BLOCK*x_bytes);

        hwpf_param_x_rem.stride = 1;
        hwpf_param_x_rem.nlines = 1;
        hwpf_param_x_rem.nblocks = __bytes_to_cachelines(PREFETCH_REMAINING*x_bytes);

        hwpf_throttle.nwait = 4;
        hwpf_throttle.ninflight = 0;

        hwpf_set_params(HWPF_ENGINE_X_ID, &hwpf_param_x);
        hwpf_set_throttle(HWPF_ENGINE_X_ID, &hwpf_throttle);

        if (PREFETCH_BLOCKS) {
            // We use REM engines to prefetch first block while we use it
            hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x);
            hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);
        } else {
            // Else it is also the last block !
            hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
            hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);
        }

        const uintptr_t _addr_x = _x;
        hwpf_wait(HWPF_ENGINE_X_REM_ID);
        hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);

        if (PREFETCH_REMAINING) {
            // If the first block is not also the last one, we prepare prefetch
            // of the last one
            hwpf_set_params(HWPF_ENGINE_X_REM_ID, &hwpf_param_x_rem);
            hwpf_set_throttle(HWPF_ENGINE_X_REM_ID, &hwpf_throttle);
        }

        for (int i = 0; i < PREFETCH_BLOCKS - 1; i++) {
            //  Prefetch next block
            const uintptr_t _addr_x = _x + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_ID, _addr_x, HWPF_ENABLE);

            int elems = PREFETCH_ELEM_PER_BLOCK;
            __vcopy_16_iter(&_x, &_y, &elems, EVP0, EVP1, x_bytes, y_bytes);
            n -= PREFETCH_ELEM_PER_BLOCK;
        }

        if (PREFETCH_REMAINING) {
            //  Prefetch last block
            const uintptr_t _addr_x = _x + PREFETCH_ELEM_PER_BLOCK*x_bytes;
            hwpf_wait(HWPF_ENGINE_X_REM_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_X_REM_ID, _addr_x, HWPF_ENABLE);
        }

        int elems = PREFETCH_ELEM_PER_BLOCK;
        __vcopy_16_iter(&_x, &_y, &elems, EVP0, EVP1, x_bytes, y_bytes);
        n -= PREFETCH_ELEM_PER_BLOCK;
    }
#endif

    __vcopy_16_iter(&_x, &_y, &n, EVP0, EVP1, x_bytes, y_bytes);

    while (n >= 8) {
        __vcopy_8(_x, _y, EVP0, EVP1);
        n  -= 8;
        _x += 8*x_bytes;
        _y += 8*y_bytes;
    }
    while (n >= 4) {
        __vcopy_4(_x, _y, EVP0, EVP1);
        n  -= 4;
        _x += 4*x_bytes;
        _y += 4*y_bytes;
    }
    while (n > 0) {
        __vcopy_1(_x, _y, EVP0, EVP1);
        n  -= 1;
        _x += 1*x_bytes;
        _y += 1*y_bytes;
    }

    /* restore environment */
    pser_evp(old_evp0, EVP0);
    pser_evp(old_evp1, EVP1);

    VBLASPERFMONITOR_FUNCTION_END;

}

void vcopy_d_v(int n, const double *x, vpfloat_off_t x_inc,
               void *y, vpfloat_evp_t y_evp, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vcopy(n, x, VPFLOAT_EVP_DOUBLE, y, y_evp, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;

}

void vcopy_v_d(int n, const void *x, vpfloat_evp_t x_evp,
               double *y, vpfloat_off_t y_inc, char enable_prefetch)
{
    VBLASPERFMONITOR_FUNCTION_BEGIN;

    vcopy(n, x, x_evp, y, VPFLOAT_EVP_DOUBLE, enable_prefetch);

    VBLASPERFMONITOR_FUNCTION_END;
}
