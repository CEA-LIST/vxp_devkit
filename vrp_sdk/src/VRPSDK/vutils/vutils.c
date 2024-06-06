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
 *  @file        vutils.c
 *  @author      Cesar Fuguet-Tortolero
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "VRPSDK/vutils.h"
#include "VRPSDK/vmath.h"

void vprint_vector(const char* nm, int n, void *v, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_efp0 = pger_efp(EFP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);
    for (int i = 0; i < n; i++) {
        ple(P0, (uintptr_t)v + i*VPFLOAT_SIZEOF(env), 0, EVP0);
        double val = rawtod(pcvt_p_d(P0, EFP0));
        printf("%s[%d] = %f\n", nm, i, val);
    }

    // restore environment
    pser_efp(old_efp0, EFP0);
    pser_evp(old_evp0, EVP0);
}

void vprint_matrix(const char* nm, int m, int n, int ld,
                   void *mat, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_efp0 = pger_efp(EFP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ple(P0, (uintptr_t)mat + (i*ld + j)*VPFLOAT_SIZEOF(env), 0, EVP0);
            double val = rawtod(pcvt_p_d(P0, EFP0));
            printf("%s[%d][%d] = %f\n", nm, i, j, val);
        }
    }

    // restore environment
    pser_efp(old_efp0, EFP0);
    pser_evp(old_evp0, EVP0);
}

double vtod(const void *v, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_efp0 = pger_efp(EFP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);
    ple(P0, (uintptr_t)v, 0, EVP0);
    double ret = rawtod(pcvt_p_d(P0, EFP0));

    // restore environment
    pser_efp(old_efp0, EFP0);
    pser_evp(old_evp0, EVP0);

    return ret;
}

void dtov(double d, void *v, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pcvt_d_p(P0, dtoraw(d));
    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pse(P0, (uintptr_t)v, 0, EVP0);

    // restore environment
    pser_evp(old_evp0, EVP0);
}

float vtof(const void *v, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);
    const uint64_t old_efp0 = pger_efp(EFP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pser_efp(pack_efp(1, vpfloat_get_rm_mem()), EFP0);
    ple(P0, (uintptr_t)v, 0, EVP0);
    float ret = rawtof(pcvt_p_f(P0, EFP0));

    // restore environment
    pser_efp(old_efp0, EFP0);
    pser_evp(old_evp0, EVP0);

    return ret;
}

void ftov(float f, void *v, vpfloat_evp_t env)
{
    // save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pcvt_f_p(P0, ftoraw(f));
    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    pse(P0, (uintptr_t)v, 0, EVP0);

    // restore environment
    pser_evp(old_evp0, EVP0);
}

void * vbuild(const uint64_t mantissa_1st_chunk, const uint64_t exponent, char a_sign, void *r, vpfloat_evp_t r_env) {
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    //  allocate memory for the result if it is not already done
    if (r == NULL) r = (void*)malloc(VPFLOAT_SIZEOF(r_env));    

    // Build the vpfloat from function arguments
    pmv_x_p_exp(P0, exponent);
    pmv_x_p_chunk(P0, 0, mantissa_1st_chunk);
    pmv_x_p_sign(P0, a_sign);

    //  initialize environments
    const vpfloat_rm_t RM = vpfloat_get_rm_mem();
    pser_evp(pack_evp(r_env.bis, RM, r_env.es, r_env.stride), EVP0);
    
    //  write-back the built number
    pse(P0, (uintptr_t)r, 0, EVP0);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return r;
}

#define EXTRACT_CHUNKS_TO_STRING(reg_name)                                      \
    case reg_name:                                                              \
                                                                                \
        l_chunk_value = pmv_p_x_chunk(0, reg_name);                             \
        sprintf(l_buffer_offset, "%016lx", l_chunk_value);                      \
        l_buffer_offset += 16;                                                  \
                                                                                \
        if ( l_chunk_count >= 1) {                                              \
            l_chunk_value = pmv_p_x_chunk(1, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 2) {                                              \
            l_chunk_value = pmv_p_x_chunk(2, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 3) {                                              \
            l_chunk_value = pmv_p_x_chunk(3, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 4) {                                              \
            l_chunk_value = pmv_p_x_chunk(4, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 5) {                                              \
            l_chunk_value = pmv_p_x_chunk(5, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 6) {                                              \
            l_chunk_value = pmv_p_x_chunk(6, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        if ( l_chunk_count >= 7) {                                              \
            l_chunk_value = pmv_p_x_chunk(7, reg_name);                         \
            sprintf(l_buffer_offset, "%016lx", l_chunk_value);                  \
            l_buffer_offset += 16;                                              \
        }                                                                       \
                                                                                \
        break;                                                                  

#define EXTRACT_TOSTRING_DATA(reg_name)                         \
    case reg_name:                                              \
        l_chunk_count = (uint16_t)pmv_p_x_len(reg_name);        \
        l_sign = pmv_p_x_sign(reg_name);                        \
        l_sum = pmv_p_x_sum(0x0000000F, reg_name);              \
        l_exp = pmv_p_x_exp(reg_name);                          \
        break;

char * vregtostring(const unsigned int preg, const unsigned int ev) {
    char * l_buffer;
    char * l_buffer_offset;
    uint16_t l_chunk_count;
    uint64_t l_chunk_value, l_sign, l_sum, l_exp;
    uint8_t l_one_flag_is_set = 0;

    // Allocate memory to store pnumber string representation
    l_buffer = (char *)malloc(sizeof(char) * 1024);

    if ( l_buffer == NULL) {
        return NULL;
    }

    memset(l_buffer, 0, 1024);

    l_buffer_offset = l_buffer;

    switch(preg) {
        EXTRACT_TOSTRING_DATA(P0)
        EXTRACT_TOSTRING_DATA(P1)
        EXTRACT_TOSTRING_DATA(P2)
        EXTRACT_TOSTRING_DATA(P3)
        EXTRACT_TOSTRING_DATA(P4)
        EXTRACT_TOSTRING_DATA(P5)
        EXTRACT_TOSTRING_DATA(P6)
        EXTRACT_TOSTRING_DATA(P7)
        EXTRACT_TOSTRING_DATA(P8)
        EXTRACT_TOSTRING_DATA(P9)
        EXTRACT_TOSTRING_DATA(P10)
        EXTRACT_TOSTRING_DATA(P11)
        EXTRACT_TOSTRING_DATA(P12)
        EXTRACT_TOSTRING_DATA(P13)
        EXTRACT_TOSTRING_DATA(P14)
        EXTRACT_TOSTRING_DATA(P15)
        EXTRACT_TOSTRING_DATA(P16)
        EXTRACT_TOSTRING_DATA(P17)
        EXTRACT_TOSTRING_DATA(P18)
        EXTRACT_TOSTRING_DATA(P19)
        EXTRACT_TOSTRING_DATA(P20)
        EXTRACT_TOSTRING_DATA(P21)
        EXTRACT_TOSTRING_DATA(P22)
        EXTRACT_TOSTRING_DATA(P23)
        EXTRACT_TOSTRING_DATA(P24)
        EXTRACT_TOSTRING_DATA(P25)
        EXTRACT_TOSTRING_DATA(P26)
        EXTRACT_TOSTRING_DATA(P27)
        EXTRACT_TOSTRING_DATA(P28)
        EXTRACT_TOSTRING_DATA(P29)
        EXTRACT_TOSTRING_DATA(P30)
        EXTRACT_TOSTRING_DATA(P31)
    }


    /***************************************
     * Display flags
     ***************************************/
    if ( l_sum != 0 ) {
        sprintf(l_buffer_offset, "Flags=");
        l_buffer_offset += 6;

        if ( l_sum & 0x1  ) {
            sprintf(l_buffer_offset, "zero");
            l_buffer_offset += 4;
            l_one_flag_is_set++;
        }

        if ( ( l_sum >> 1 ) & 0x1 ) {
            if ( l_one_flag_is_set ) {
                sprintf(l_buffer_offset, ",");    
                l_buffer_offset += 1;
            }
            sprintf(l_buffer_offset, "inf");
            l_buffer_offset += 3;
            l_one_flag_is_set++;
        }

        if ( ( l_sum >> 2 ) & 0x1 ) {
            if ( l_one_flag_is_set ) {
                sprintf(l_buffer_offset, ",");    
                l_buffer_offset += 1;
            }        
            sprintf(l_buffer_offset, "qNaN");
            l_buffer_offset += 4;    
            l_one_flag_is_set++;
        }

        if ( ( l_sum >> 3 ) & 0x1 ) {
            if ( l_one_flag_is_set ) {
                sprintf(l_buffer_offset, ",");    
                l_buffer_offset += 1;
            }        
            sprintf(l_buffer_offset, "sNaN");
            l_buffer_offset += 4;    
            l_one_flag_is_set++;
        }

        sprintf(l_buffer_offset, " ");
        l_buffer_offset += 1;
    }

    /***************************************
     * Display data
     ***************************************/
    if ( l_sum == 0 ) {
        sprintf(l_buffer_offset, "E=0x%lx ", l_exp);
        l_buffer_offset = l_buffer + strlen(l_buffer);

        sprintf(l_buffer_offset, "L=%d ", l_chunk_count);
        l_buffer_offset = l_buffer + strlen(l_buffer);

        sprintf(l_buffer_offset, "M=%c", l_sign == 0 ? '+' : '-');
        l_buffer_offset = l_buffer + strlen(l_buffer);

        switch(preg) {
            EXTRACT_CHUNKS_TO_STRING(P0)
            EXTRACT_CHUNKS_TO_STRING(P1)
            EXTRACT_CHUNKS_TO_STRING(P2)
            EXTRACT_CHUNKS_TO_STRING(P3)
            EXTRACT_CHUNKS_TO_STRING(P4)
            EXTRACT_CHUNKS_TO_STRING(P5)
            EXTRACT_CHUNKS_TO_STRING(P6)
            EXTRACT_CHUNKS_TO_STRING(P7)
            EXTRACT_CHUNKS_TO_STRING(P8)
            EXTRACT_CHUNKS_TO_STRING(P9)
            EXTRACT_CHUNKS_TO_STRING(P10)
            EXTRACT_CHUNKS_TO_STRING(P11)
            EXTRACT_CHUNKS_TO_STRING(P12)
            EXTRACT_CHUNKS_TO_STRING(P13)
            EXTRACT_CHUNKS_TO_STRING(P14)
            EXTRACT_CHUNKS_TO_STRING(P15)
            EXTRACT_CHUNKS_TO_STRING(P16)
            EXTRACT_CHUNKS_TO_STRING(P17)
            EXTRACT_CHUNKS_TO_STRING(P18)
            EXTRACT_CHUNKS_TO_STRING(P19)
            EXTRACT_CHUNKS_TO_STRING(P20)
            EXTRACT_CHUNKS_TO_STRING(P21)
            EXTRACT_CHUNKS_TO_STRING(P22)
            EXTRACT_CHUNKS_TO_STRING(P23)
            EXTRACT_CHUNKS_TO_STRING(P24)
            EXTRACT_CHUNKS_TO_STRING(P25)
            EXTRACT_CHUNKS_TO_STRING(P26)
            EXTRACT_CHUNKS_TO_STRING(P27)
            EXTRACT_CHUNKS_TO_STRING(P28)
            EXTRACT_CHUNKS_TO_STRING(P29)
            EXTRACT_CHUNKS_TO_STRING(P30)
            EXTRACT_CHUNKS_TO_STRING(P31)
        }
    }

    *(l_buffer_offset) = '\0';

    return l_buffer;
}

char * vtostring(const void *v, vpfloat_evp_t env) {
    char * l_buffer;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    ple(P0, (uintptr_t)v, 0, EVP0);

    l_buffer = vregtostring(P0, EVP0);

    char * l_buffer_full = (char *)malloc(strlen(l_buffer) + 20 );

    /* FIXME (cf247698): this is wrong. Behaviour is undefined if the
     * destination is also one of the sources */
    sprintf(l_buffer_full, "%s - %p", l_buffer, v);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    free(l_buffer);

    return l_buffer_full;
}

/**
 * vblas_vtrans
 * by Valentin Isaac--Chassande
 *
 * Matrix transposition without copy
 *
 * Use bitset.h (common/bitset.h)
 */
static void vtrans_square_1(int n,
                            const void *a, int a_bytes)
{
    uintptr_t a_ptr;
    uintptr_t at_ptr;
    int a_row_bytes = n * a_bytes;

    if  (n > 1) {                                    /* if n=1, there is nothing to do */
        for (int i = 1; i < n; i++) {
            a_ptr = (uintptr_t)a + i * a_row_bytes;  /* select a line */
            at_ptr = a_ptr - i*(n-1);                /* select the corresponding column */
            for (int j = 0; j < i; j++) {
                ple(P0, a_ptr, 0, EVP0);             /* operate the exchange */
                ple(P1, at_ptr, 0, EVP0);
                pse(P0, at_ptr, 0, EVP0);
                pse(P1, a_ptr, 0, EVP0);
                a_ptr += a_bytes;                    /* move onto the line */
                at_ptr += a_row_bytes;               /* move onto the column */
            }
        }
    }
}

static void vtrans_rect_1(int n, int p,
                          const void *a, int a_bytes)
{
    uintptr_t a_ptr;
    uintptr_t at_ptr;
    int a_row_bytes = p * a_bytes;
    int pos, new_i, new_j;
    bitset_t * bitset = (bitset_t *) malloc(sizeof(bitset_t));

    if  (n > 1 && p > 1) {                           /* if n=1 or p=1, there is nothing to do */
        bitset_init(bitset, n*p);                    /* allocation and initialization of a */
                                                     /*     matrice of bits                */
        for (int i = 1; i < n; i++) {
            a_ptr = (uintptr_t)a + i * a_row_bytes;  /* select a line */
            for (int j = 0; j < p; j++) {
                pos = i*p+j;
                if (!bitset_test(bitset, pos)) {     /* check if the element has already */
                                                     /*     been compute                 */
                    new_i = pos%n;
                    new_j = pos/n;
                    at_ptr = (uintptr_t)a + new_i*p+new_j; /* select the exchange position */
                    if  (a_ptr != at_ptr) {          /* check if the two positions are the same */
                        ple(P0, a_ptr, 0, EVP0);     /* operate the exchange */
                        ple(P1, at_ptr, 0, EVP0);
                        pse(P0, at_ptr, 0, EVP0);
                        pse(P1, a_ptr, 0, EVP0);
                        bitset_set(bitset, new_i*p+new_j); /* mark the exchange position */
                    }
                }
                a_ptr += a_bytes;                    /* move onto the line */
            }
        }
        bitset_destroy(bitset);                      /* free the matrice of bits */
    }

    free((void *) bitset);
}

static void vtrans_square_4(int n,
                            const void *a, int a_bytes)
{
    uintptr_t a_ptr;
    uintptr_t at_ptr;
    int a_row_bytes = n * a_bytes;

    if  (n > 1) {
        for (int i = 1; i < n; i++) {
            a_ptr = (uintptr_t)a + i * a_row_bytes;
            at_ptr = a_ptr - i*(n-1);
            int j = i;
            while (j >= 4) {                         /* four exchanges at once */
                ple(P0, a_ptr+0*a_bytes, 0, EVP0);
                ple(P2, a_ptr+1*a_bytes, 0, EVP0);
                ple(P4, a_ptr+2*a_bytes, 0, EVP0);
                ple(P6, a_ptr+3*a_bytes, 0, EVP0);

                ple(P1, at_ptr+0*a_row_bytes, 0, EVP0);
                ple(P3, at_ptr+1*a_row_bytes, 0, EVP0);
                ple(P5, at_ptr+2*a_row_bytes, 0, EVP0);
                ple(P7, at_ptr+3*a_row_bytes, 0, EVP0);

                pse(P0, at_ptr+0*a_row_bytes, 0, EVP0);
                pse(P2, at_ptr+1*a_row_bytes, 0, EVP0);
                pse(P4, at_ptr+2*a_row_bytes, 0, EVP0);
                pse(P6, at_ptr+3*a_row_bytes, 0, EVP0);

                pse(P1, a_ptr+0*a_bytes, 0, EVP0);
                pse(P3, a_ptr+1*a_bytes, 0, EVP0);
                pse(P5, a_ptr+2*a_bytes, 0, EVP0);
                pse(P7, a_ptr+3*a_bytes, 0, EVP0);

                a_ptr += 4*a_bytes;
                at_ptr += 4*a_row_bytes;
                j-=4;
            }
            if (j >= 2) {                            /* two exchanges at once */
                ple(P0, a_ptr+0*a_bytes, 0, EVP0);
                ple(P2, a_ptr+1*a_bytes, 0, EVP0);

                ple(P1, at_ptr+0*a_row_bytes, 0, EVP0);
                ple(P3, at_ptr+1*a_row_bytes, 0, EVP0);

                pse(P0, at_ptr+0*a_row_bytes, 0, EVP0);
                pse(P2, at_ptr+1*a_row_bytes, 0, EVP0);

                pse(P1, a_ptr+0*a_bytes, 0, EVP0);
                pse(P3, a_ptr+1*a_bytes, 0, EVP0);

                a_ptr += 2*a_bytes;
                at_ptr += 2*a_row_bytes;
                j-=2;
            }
            if (j == 1) {                            /* one exchange */
                ple(P0, a_ptr, 0, EVP0);
                ple(P1, at_ptr, 0, EVP0);

                pse(P0, at_ptr, 0, EVP0);
                pse(P1, a_ptr, 0, EVP0);
            }
        }
    }
}

void vtrans(int precision, int n, int p,
            const void *a, vpfloat_evp_t a_evp)
{
    /* save environment */
    const uint64_t old_evp0 = pger_evp(EVP0);

    /* set compute and memory environments */
    pser_evp(pack_evp(a_evp.bis, vpfloat_get_rm_mem(),
                      a_evp.es, a_evp.stride),
             EVP0);

    /* execute transpose */
    if (n != p) {
        vtrans_rect_1(n, p, a, VPFLOAT_SIZEOF(a_evp));
    } else {
        vtrans_square_4(n, a, VPFLOAT_SIZEOF(a_evp));
    }

    if (0) vtrans_square_1(n, a, VPFLOAT_SIZEOF(a_evp));

    /* restore environment */
    pser_evp(old_evp0, EVP0);
}

uint64_t vlen(const void *v, vpfloat_evp_t env) {
    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    ple(P0, (uintptr_t)v, 0, EVP0);

    uint64_t l_len = pmv_p_x_len(P0);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return l_len;
}

uint64_t isNaN(const void *v, vpfloat_evp_t env) {
    uint64_t in_n_isNaNquiet_int , in_n_isNaNsignaling_int;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    ple(P0, (uintptr_t)v, 0, EVP0);

    in_n_isNaNquiet_int     = pmv_p_x_sum(SUMM_IS_QNAN, P0);
    in_n_isNaNsignaling_int = pmv_p_x_sum(SUMM_IS_SNAN, P0);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return (in_n_isNaNquiet_int != 0 || in_n_isNaNsignaling_int != 0);
}

uint64_t isInf(const void *v, vpfloat_evp_t env) {
    uint64_t in_n_isINF_int;

    //  save environment
    const uint64_t old_evp0 = pger_evp(EVP0);

    pser_evp(pack_evp(env.bis, vpfloat_get_rm_mem(), env.es, env.stride), EVP0);
    ple(P0, (uintptr_t)v, 0, EVP0);

    in_n_isINF_int = pmv_p_x_sum(SUMM_IS_INF,  P0);

    //  restore environment
    pser_evp(old_evp0, EVP0);

    return ( in_n_isINF_int != 0 );
}
