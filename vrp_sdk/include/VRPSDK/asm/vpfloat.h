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
 *  @file        vpfloat.h
 *  @author      Cesar Fuguet-Tortolero
 *  @date        October, 2021
 *  @brief       Definition of wrappers for the assembly instructions of the
 *               Variable Precision Coprocessor (RISC-V Xvpfloat instruction
 *               set extension).
 */
#ifndef __ASM_VPFLOAT_H__
#define __ASM_VPFLOAT_H__

#include <stdint.h>
#include <stddef.h>
#include "VRPSDK/asm/vpfloat_defs.h"

/*
 * Environment definitions
 */
/* Stride */
typedef int32_t        vpfloat_off_t;
typedef vpfloat_off_t  vpfloat_stride_t;

/* Precision */
typedef uint16_t       vpfloat_prec_t;
typedef vpfloat_prec_t vpfloat_bis_t;
typedef vpfloat_prec_t vpfloat_wp_t;

/* Exponent */
typedef uint16_t       vpfloat_es_t;

/* Rounding mode */
typedef uint8_t        vpfloat_rm_t;

typedef struct {
    vpfloat_stride_t stride;
    vpfloat_bis_t    bis;
    vpfloat_es_t     es;
} vpfloat_evp_t;

typedef struct {
    vpfloat_stride_t stride;
} vpfloat_efp_t;

typedef struct {
    vpfloat_wp_t     wp;
} vpfloat_ec_t;

/* Predefined environments */
static const vpfloat_evp_t VPFLOAT_EVP_DOUBLE = {
    .stride = 1,
    .bis    = 64,
    .es     = 11
};

static const vpfloat_evp_t VPFLOAT_EVP_FLOAT = {
    .stride = 1,
    .bis    = 32,
    .es     = 8
};

static const vpfloat_evp_t VPFLOAT_EVP_HALF = {
    .stride = 1,
    .bis    = 16,
    .es     = 5
};

static const vpfloat_evp_t VPFLOAT_EVP_MAX = {
    .stride = 1,
    .bis    = BISMAX,
    .es     = ELEN
};

#define VPFLOAT_IS_DOUBLE(x) (\
        (x.bis == VPFLOAT_EVP_DOUBLE.bis) && \
        (x.es  == VPFLOAT_EVP_DOUBLE.es))

#define VPFLOAT_IS_FLOAT(x) (\
        (x.bis == VPFLOAT_EVP_FLOAT.bis) && \
        (x.es  == VPFLOAT_EVP_FLOAT.es))

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Pack/unpack functions for environment registers
 */
/* Compute environment */
volatile static inline uint64_t pack_ec(vpfloat_wp_t wp,
                                        vpfloat_rm_t rm)
{
    return VPFLOAT_BITSET(RND_BITS, RND_OFFSET, (uint64_t)(rm    )) |
           VPFLOAT_BITSET( WP_BITS,  WP_OFFSET, (uint64_t)(wp - 1));
}

volatile static inline void unpack_ec(vpfloat_wp_t *wp,
                                      vpfloat_rm_t *rm,
                                      const uint64_t ec)
{
    if (wp != NULL) *wp = VPFLOAT_BITGET( WP_BITS,  WP_OFFSET, (uint64_t)ec) + 1;
    if (rm != NULL) *rm = VPFLOAT_BITGET(RND_BITS, RND_OFFSET, (uint64_t)ec);
}

/* Variable memory environment */
volatile static inline uint64_t pack_evp(vpfloat_bis_t    bis,
                                         vpfloat_rm_t     rm,
                                         vpfloat_es_t     es,
                                         vpfloat_stride_t stride)
{
    return VPFLOAT_BITSET(   BIS_BITS,    BIS_OFFSET, (uint64_t)(bis    - 1)) |
           VPFLOAT_BITSET(   RND_BITS,    RND_OFFSET, (uint64_t)(rm        )) |
           VPFLOAT_BITSET(    ES_BITS,     ES_OFFSET, (uint64_t)(es     - 1)) |
           VPFLOAT_BITSET(STRIDE_BITS, STRIDE_OFFSET, (uint64_t)(stride - 1));
}

volatile static inline void unpack_evp(vpfloat_bis_t    *bis,
                                       vpfloat_rm_t     *rm,
                                       vpfloat_es_t     *es,
                                       vpfloat_stride_t *stride,
                                       const uint64_t    evp)
{
    if (   bis != NULL) *bis    = VPFLOAT_BITGET(   BIS_BITS,    BIS_OFFSET, evp) + 1;
    if (    rm != NULL) *rm     = VPFLOAT_BITGET(   RND_BITS,    RND_OFFSET, evp);
    if (    es != NULL) *es     = VPFLOAT_BITGET(    ES_BITS,     ES_OFFSET, evp) + 1;
    if (stride != NULL) *stride = VPFLOAT_BITGET(STRIDE_BITS, STRIDE_OFFSET, evp) + 1;
}

/* Fixed memory environment */
volatile static inline uint64_t pack_efp(vpfloat_stride_t stride,
                                         vpfloat_rm_t     rm)
{
    return VPFLOAT_BITSET(       RND_BITS,        RND_OFFSET, (uint64_t)(rm        )) |
           VPFLOAT_BITSET(STRIDE_EFP_BITS, STRIDE_EFP_OFFSET, (uint64_t)(stride - 1));
}

volatile static inline void unpack_efp(vpfloat_stride_t *stride,
                                       vpfloat_rm_t     *rm,
                                       const uint64_t    efp)
{
    if (stride != NULL) *stride = VPFLOAT_BITGET(STRIDE_EFP_BITS, STRIDE_EFP_OFFSET, efp) + 1;
    if (    rm != NULL) *rm     = VPFLOAT_BITGET(       RND_BITS,        RND_OFFSET, efp);
}


/*
 * Compare environment variables
 */
static inline int vpfloat_evp_is_equal(const vpfloat_evp_t a, const vpfloat_evp_t b)
{
    return (a.stride == b.stride) &&
           (a.bis    == b.bis) &&
           (a.es     == b.es);

}

static inline int vpfloat_evp_is_double(const vpfloat_evp_t a)
{
    return vpfloat_evp_is_equal(a, VPFLOAT_EVP_DOUBLE);
}

static inline int vpfloat_evp_is_single(const vpfloat_evp_t a)
{
    return vpfloat_evp_is_equal(a, VPFLOAT_EVP_FLOAT);
}


/*
 *  Set environment register
 */
volatile static inline void pser_evp(const uint64_t val, const unsigned int env)
{
    asm volatile ("pser evp%c[env], %[val]"
                 : /* no outputs */
                 : [env]"I"(env), [val]"r"(val)
                 : "memory");
}

volatile static inline void pser_efp(const uint64_t val, const unsigned int env)
{
    asm volatile ("pser efp%c[env], %[val]"
                 : /* no outputs */
                 : [env]"I"(env), [val]"r"(val)
                 : "memory");
}

volatile static inline void pser_ec(const uint64_t val, const unsigned int env)
{
    asm volatile ("pser ec%c[env], %[val]"
                 : /* no outputs */
                 : [env]"I"(env), [val]"r"(val)
                 : "memory");
}

/*
 *  Get environment register
 */
volatile static inline uint64_t pger_evp(const unsigned int env)
{
    register uint64_t ret;

    asm volatile ("pger %[ret], evp%c[env]"
                 : [ret]"=&r"(ret)
                 : [env]"I"(env));

    return ret;
}

volatile static inline uint64_t pger_efp(const unsigned int env)
{
    register uint64_t ret;

    asm volatile ("pger %[ret], efp%c[env]"
                 : [ret]"=&r"(ret)
                 : [env]"I"(env));

    return ret;
}

volatile static inline uint64_t pger_ec(const unsigned int env)
{
    register uint64_t ret;

    asm volatile ("pger %[ret], ec%c[env]"
                 : [ret]"=&r"(ret)
                 : [env]"I"(env));

    return ret;
}

/******************************************************
 *  Conversion instructions
 *****************************************************/
#define PCVT_X_P(X, _dst, _src) \
    asm volatile ("pcvt." #X ".p p%c[dst], %[src]  \n" \
                 :                                     \
                 : [src]"r"(_src),  [dst]"I"(_dst))

#define PCVT_P_X(X, _src, _ret, _env) \
    asm volatile ("pcvt.p." #X " %[dst], p%c[src], efp%c[env]  \n" \
                 : [dst]"=&r"(_ret)                                \
                 : [src]"I"(_src), [env]"I"(_env))

volatile static inline void pcvt_h_p(const unsigned dst, const uint16_t src)
{
    PCVT_X_P(h, dst, src);
}

volatile static inline void pcvt_f_p(const unsigned dst, const uint32_t src)
{
    PCVT_X_P(f, dst, src);
}

volatile static inline void pcvt_d_p(const unsigned dst, const uint64_t src)
{
    PCVT_X_P(d, dst, src);
}

volatile static inline uint16_t pcvt_p_h(const unsigned src, const unsigned env)
{
    uint16_t ret = 0;
    PCVT_P_X(h, src, ret, env);
    return ret;
}

volatile static inline uint32_t pcvt_p_f(const unsigned src, const unsigned env)
{
    uint32_t ret = 0;
    PCVT_P_X(f, src, ret, env);
    return ret;
}

volatile static inline uint64_t pcvt_p_d(const unsigned src, const unsigned env)
{
    uint64_t ret = 0;
    PCVT_P_X(d, src, ret, env);
    return ret;
}

#undef PCVT_X_P
#undef PCVT_P_X

/******************************************************
 *  Move instructions
 *****************************************************/
volatile static inline uint64_t pack_header(uint64_t sign,
                                            uint64_t is_nan_quiet,
                                            uint64_t is_nan_signaling,
                                            uint64_t is_inf,
                                            uint64_t is_zero,
                                            uint64_t l,
                                            uint64_t exponent)
{
    return VPFLOAT_BITSET(            SIGN_BITS,             SIGN_OFFSET,             sign) |
           VPFLOAT_BITSET(    IS_NAN_QUIET_BITS,     IS_NAN_QUIET_OFFSET,     is_nan_quiet) |
           VPFLOAT_BITSET(IS_NAN_SIGNALING_BITS, IS_NAN_SIGNALING_OFFSET, is_nan_signaling) |
           VPFLOAT_BITSET(          IS_INF_BITS,           IS_INF_OFFSET,           is_inf) |
           VPFLOAT_BITSET(         IS_ZERO_BITS,          IS_ZERO_OFFSET,          is_zero) |
           VPFLOAT_BITSET(               L_BITS,                L_OFFSET,                l) |
           VPFLOAT_BITSET(             EXP_BITS,              EXP_OFFSET,         exponent);
}

#define PMV_P_X(X, src) ({                         \
    uint64_t ret;                                  \
    asm volatile ("pmv." #X ".x %[dst], p%c[src]"  \
                 : [dst]"=&r"(ret)                 \
                 : [src]"I"(src));                 \
    ret;                                           \
    })

#define PMV_X_P(X, _dst, _src) \
    asm volatile ("pmv.x." #X " p%c[dst], %[src]"  \
                 :                                 \
                 : [src]"r"(_src),  [dst]"I"(_dst))

volatile static inline uint64_t pmv_p_x_chunk(const uint16_t chunk, const unsigned src)
{
    uint64_t ret;
    asm volatile ("pmv.p.x %[dst], p%c[src], %c[chunk]"
                 : [dst]"=&r"(ret)
                 : [src]"I"(src), [chunk]"I"(chunk));
    return ret;
}

volatile static inline uint64_t pmv_p_x_sum(const uint16_t mask, const unsigned src)
{
    uint64_t ret;
    asm volatile ("pmv.psum.x %[dst], p%c[src], %c[mask]"
                 : [dst]"=&r"(ret)
                 : [src]"I"(src), [mask]"I"(mask));
    return ret;
}

volatile static inline uint64_t pmv_p_x_sign(const unsigned src)
{
    return PMV_P_X(psgn, src);
}

volatile static inline uint64_t pmv_p_x_len(const unsigned src)
{
    return PMV_P_X(plen, src);
}

volatile static inline uint64_t pmv_p_x_exp(const unsigned src)
{
    return PMV_P_X(pexp, src);
}

volatile static inline void pmv_p_p(const unsigned dst, const unsigned src)
{
    asm volatile ("pmv.p.p p%c[dst], p%c[src]"
                 :
                 : [src]"I"(src), [dst]"I"(dst));
}

volatile static inline void pmv_x_p_chunk(const unsigned dst, const uint16_t chunk, const uint64_t src)
{
    asm volatile ("pmv.x.p p%c[dst], %[src], %c[chunk]"
                 :
                 : [dst]"I"(dst), [src]"r"(src), [chunk]"I"(chunk));
}

volatile static inline void pmv_x_p_sign(const unsigned dst, const uint64_t src)
{
    PMV_X_P(psgn, dst, src);
}

volatile static inline void pmv_x_p_sum(const unsigned dst, const uint64_t src)
{
    PMV_X_P(psum, dst, src);
}

volatile static inline void pmv_x_p_len(const unsigned dst, const uint64_t src)
{
    PMV_X_P(plen, dst, src);
}

volatile static inline void pmv_x_p_exp(const unsigned dst, const uint64_t src)
{
    PMV_X_P(pexp, dst, src);
}

#undef PMV_X_P
#undef PMV_P_X

/******************************************************
 *  Comparison instructions
 *****************************************************/
#define PCMP(X, _rs1, _rs2) ({                            \
    uint64_t ret;                                         \
    asm volatile("pcmp." #X " %[dst], p%c[rs1], p%c[rs2]" \
        : [dst]"=&r"(ret)                                 \
        : [rs1]"I"(_rs1), [rs2]"I"(_rs2));                \
    ret;                                                  \
    })

volatile static inline int pcmp(const unsigned rs1, const unsigned rs2)
{
    uint64_t ret;
    asm volatile("pcmp %[dst], p%c[rs1], p%c[rs2]"
                : [dst]"=&r"(ret)
                : [rs1]"I"(rs1), [rs2]"I"(rs2));
    return ret;
}

volatile static inline int pcmp_eq(const unsigned rs1, const unsigned rs2)
{
    return PCMP(eq, rs1, rs2);
}

volatile static inline int pcmp_neq(const unsigned rs1, const unsigned rs2)
{
    return PCMP(neq, rs1, rs2);
}

volatile static inline int pcmp_gt(const unsigned rs1, const unsigned rs2)
{
    return PCMP(gt, rs1, rs2);
}

volatile static inline int pcmp_lt(const unsigned rs1, const unsigned rs2)
{
    return PCMP(lt, rs1, rs2);
}

volatile static inline int pcmp_geq(const unsigned rs1, const unsigned rs2)
{
    return PCMP(geq, rs1, rs2);
}

volatile static inline int pcmp_leq(const unsigned rs1, const unsigned rs2)
{
    return PCMP(leq, rs1, rs2);
}

#undef PCMP

/******************************************************
 *  Memory instructions
 *****************************************************/
/*
 *  Load operations
 */
#define PLX(X, _rd, _addr, _idx, _env) \
    asm volatile (                                                            \
            "pl" #X " p%c[rd], efp%c[env], %[addr], %c[idx]"                  \
            : /* no outputs */                                                \
            :  [rd]"I"(_rd), [env]"I"(_env), [addr]"r"(_addr), [idx]"I"(_idx) \
            : "memory")

volatile static inline void plh(const unsigned int rd,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PLX(h, rd, addr, idx, env);
}

volatile static inline void plw(const unsigned int rd,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PLX(w, rd, addr, idx, env);
}

volatile static inline void pld(const unsigned int rd,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PLX(d, rd, addr, idx, env);
}

volatile static inline void ple(const unsigned int rd,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    asm volatile ("ple p%c[rd], evp%c[env], %[addr], %c[idx]"
                 : /* no outputs */
                 :  [rd]"I"(rd), [env]"I"(env), [addr]"r"(addr), [idx]"I"(idx)
                 : "memory");
}

#undef PLX

/*
 *  Store operations
 */
#define PSX(X, _rs1, _addr, _idx, _env)                         \
    asm volatile (                                              \
            "ps" #X " p%c[rs1], efp%c[env], %[addr], %c[idx]"   \
            : /* no outputs */                                  \
            : [rs1]"I"(_rs1), [env]"I"(_env), [addr]"r"(_addr), \
              [idx]"I"(_idx)                                    \
            : "memory");

volatile static inline void psh(const unsigned int rs1,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PSX(h, rs1, addr, idx, env);
}

volatile static inline void psw(const unsigned int rs1,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PSX(w, rs1, addr, idx, env);
}

volatile static inline void psd(const unsigned int rs1,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    PSX(d, rs1, addr, idx, env);
}

volatile static inline void pse(const unsigned int src,
                                const uintptr_t    addr,
                                const unsigned int idx,
                                const unsigned int env)
{
    asm volatile ("pse p%c[src], evp%c[env], %[addr], %c[idx]"
                 : /* no outputs */
                 : [src]"I"(src), [env]"I"(env), [addr]"r"(addr),
                   [idx]"I"(idx)
                 : "memory");
}

#undef PSX

/******************************************************
 *  Arithmetic instructions
 *****************************************************/
volatile static inline void padd(const unsigned int rd,
                                 const unsigned int rs1,
                                 const unsigned int rs2,
                                 const unsigned int ec)
{
    asm volatile("padd p%c[rd], p%c[rs1], p%c[rs2], ec%c[ec]"
                : /* no outputs */
                : [rd]"I"(rd), [rs1]"I"(rs1), [rs2]"I"(rs2),
                  [ec]"I"(ec));
}

volatile static inline void psub(const unsigned int rd,
                                 const unsigned int rs1,
                                 const unsigned int rs2,
                                 const unsigned int ec)
{
    asm volatile("psub p%c[rd], p%c[rs1], p%c[rs2], ec%c[ec]"
                : /* no outputs */
                : [rd]"I"(rd), [rs1]"I"(rs1), [rs2]"I"(rs2),
                  [ec]"I"(ec));
}

volatile static inline void pmul(const unsigned int rd,
                                 const unsigned int rs1,
                                 const unsigned int rs2,
                                 const unsigned int ec)
{
    asm volatile("pmul p%c[rd], p%c[rs1], p%c[rs2], ec%c[ec]"
                : /* no outputs */
                : [rd]"I"(rd), [rs1]"I"(rs1), [rs2]"I"(rs2),
                  [ec]"I"(ec));
}

volatile static inline void prnd(const unsigned int rd,
                                 const unsigned int rs1,
                                 const unsigned int ec)
{
    asm volatile("prnd p%c[rd], p%c[rs1], ec%c[ec]"
                : /* no outputs */
                : [rd]"I"(rd), [rs1]"I"(rs1),
                  [ec]"I"(ec));
}

#ifdef __cplusplus
}
#endif

#endif /* __ASM_VPFLOAT_H__ */
