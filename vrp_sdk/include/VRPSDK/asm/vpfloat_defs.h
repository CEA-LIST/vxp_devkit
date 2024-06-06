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
 *  @file        vpfloat_defs.h
 *  @author      Cesar Fuguet-Tortolero
 *  @date        October, 2021
 *  @brief       Definition of constants and macros for the Variable Precision
 *               Coprocessor (RISC-V Xvpfloat instruction set extension).
 */
#ifndef __ASM_VPFLOAT_DEFS_H__
#define __ASM_VPFLOAT_DEFS_H__

/* Variable-precision registers' index */
#define P0            0
#define P1            1
#define P2            2
#define P3            3
#define P4            4
#define P5            5
#define P6            6
#define P7            7
#define P8            8
#define P9            9
#define P10           10
#define P11           11
#define P12           12
#define P13           13
#define P14           14
#define P15           15
#define P16           16
#define P17           17
#define P18           18
#define P19           19
#define P20           20
#define P21           21
#define P22           22
#define P23           23
#define P24           24
#define P25           25
#define P26           26
#define P27           27
#define P28           28
#define P29           29
#define P30           30
#define P31           31

/* Environment definitions */
#define EVP0          0
#define EVP1          1
#define EVP2          2
#define EVP3          3
#define EVP4          4
#define EVP5          5
#define EVP6          6
#define EVP7          7

#define EFP0          0
#define EFP1          1
#define EFP2          2
#define EFP3          3
#define EFP4          4
#define EFP5          5
#define EFP6          6
#define EFP7          7

#define EC0           0
#define EC1           1
#define EC2           2
#define EC3           3
#define EC4           4
#define EC5           5
#define EC6           6
#define EC7           7

/* Rounding modes */
#define RNE           0
#define RTZ           1
#define RDN           2
#define RUP           3
#define RMM           4

/* Chunks identifier */
#define HEADER_CHUNK  31
#define P_CHUNK0      0
#define P_CHUNK1      1
#define P_CHUNK2      2
#define P_CHUNK3      3
#define P_CHUNK4      4
#define P_CHUNK5      5
#define P_CHUNK6      6
#define P_CHUNK7      7
#define P_CHUNK_MAX_COUNT   8

#define P_CHUNK_LEN   64

/* Utility macros */
#define BYS(bis)                       (((bis) + 7)/8)
#define VPFLOAT_MASK(bits)             ((1ULL << (bits)) - 1)
#define VPFLOAT_BITSET(bits, off, val) (((val) & VPFLOAT_MASK(bits)) << (off))
#define VPFLOAT_BITGET(bits, off, val) (((val) >> (off)) & VPFLOAT_MASK(bits))
#define VPFLOAT_SIZEOF(env)            (((env).stride)*BYS((env).bis))

// Maximum values
#define ELEN   18  // internal exponent length
#define MLEN   512 // internal mantissa length
#define LLEN   3   // length of the L header field, it contains the number of valid 64-bit mantissa chunks (-1) of the VP number
#define HLEN   26  // P register header size (bits)
#define VPLEN  538 // Total size of internal scalar P register (bits)
#define BISMAX 531 // Maximum size in memory (bits)
#define BYSMAX 67  // Maximum size in memory (bytes)
#define WP_MAX (MLEN-1)

#define IL_ESM1_MIN 0
#define IL_ESM1_MAX 17

/* Header configuration */
#define SIGN_BITS                 1
#define SIGN_OFFSET               25
#define SIGN_MASK                 (VPFLOAT_MASK(SIGN_BITS) << (SIGN_OFFSET))
#define IS_NAN_SIGNALING_BITS     1
#define IS_NAN_SIGNALING_OFFSET   24
#define IS_NAN_SIGNALING_MASK     (VPFLOAT_MASK(IS_NAN_SIGNALING_BITS) << (IS_NAN_SIGNALING_OFFSET))
#define IS_NAN_QUIET_BITS         1
#define IS_NAN_QUIET_OFFSET       23
#define IS_NAN_QUIET_MASK         (VPFLOAT_MASK(IS_NAN_QUIET_BITS) << (IS_NAN_QUIET_OFFSET))
#define IS_INF_BITS               1
#define IS_INF_OFFSET             22
#define IS_INF_MASK               (VPFLOAT_MASK(IS_INF_BITS) << (IS_INF_OFFSET))
#define IS_ZERO_BITS              1
#define IS_ZERO_OFFSET            21
#define IS_ZERO_MASK              (VPFLOAT_MASK(IS_ZERO_BITS) << (IS_ZERO_OFFSET))
#define L_BITS                    3
#define L_OFFSET                  18
#define L_MASK                    (VPFLOAT_MASK(L_BITS) << (L_OFFSET))
#define EXP_BITS                  18
#define EXP_OFFSET                0
#define EXP_MASK                  (VPFLOAT_MASK(EXP_BITS) << (EXP_OFFSET))

/* Summary bits masks (for usage with the mv.sum instruction) */
#define SUMM_IS_SNAN              0x8
#define SUMM_IS_QNAN              0x4
#define SUMM_IS_INF               0x2
#define SUMM_IS_ZERO              0x1

/* Environment configuration */
#define WP_BITS           16
#define WP_OFFSET         0
#define WP_MASK           (VPFLOAT_MASK(WP_BITS) << (WP_OFFSET))
#define BIS_BITS          16
#define BIS_OFFSET        0
#define BIS_MASK          (VPFLOAT_MASK(BIS_BITS) << (BIS_OFFSET))
#define RND_BITS          3
#define RND_OFFSET        16
#define RND_MASK          (VPFLOAT_MASK(RND_BITS) << (RND_OFFSET))
#define ES_BITS           8
#define ES_OFFSET         24
#define ES_MASK           (VPFLOAT_MASK(ES_BITS) << (ES_OFFSET))
#define STRIDE_BITS       16
#define STRIDE_OFFSET     32
#define STRIDE_MASK       (VPFLOAT_MASK(STRIDE_BITS) << (STRIDE_OFFSET))
#define STRIDE_EFP_BITS   16
#define STRIDE_EFP_OFFSET 0
#define STRIDE_EFP_MASK   (VPFLOAT_MASK(STRIDE_EFP_BITS) << (STRIDE_EFP_OFFSET))

#endif /* __ASM_VPFLOAT_DEFS_H__ */
