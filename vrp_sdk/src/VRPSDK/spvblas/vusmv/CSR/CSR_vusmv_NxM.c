/**
* Copyright 2023 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        CSR_vusmv_NxM.c
 *  @author      Jerome Fereyre
 */

#include <stdint.h>
#include <stdio.h>
#include <VRPSDK/spvblas.h>
#include "VRPSDK/spvblas/vusmv/CSR_vusmv_NxM.h"
#include "VRPSDK/asm/vpfloat.h"
#include "VRPSDK/vutils.h"
#include "VRPSDK/vblas_perfmonitor.h"

#ifndef SPVBLAS_ENABLE_HWPF
#define SPVBLAS_ENABLE_HWPF 1
#endif
#if SPVBLAS_ENABLE_HWPF
#include "drivers/hwpf_dcache.h"
#endif

#define ROW_ACCU_REG 	P31
#define BETA_REG 	    P30
#define ALPHA_REG 	    P29
#define X_REG 		    P28
#define Y_REG 		    P27
#define A_REG 		    P26

#define X0_REG          P0
#define X1_REG          P1
#define X2_REG          P2
#define X3_REG          P3
#define X4_REG          P4
#define X5_REG          P5
#define X6_REG          P6
#define X7_REG          P7
#define ACC0_REG          P8
#define ACC1_REG          P9
#define ACC2_REG         P10
#define ACC3_REG         P11
#define ACC4_REG         P12
#define ACC5_REG         P13
#define ACC6_REG         P14
#define ACC7_REG         P15
#define A0_REG          P16
#define A1_REG          P17
#define A2_REG          P18
#define A3_REG          P19
#define A4_REG          P20
#define A5_REG          P21
#define A6_REG          P22
#define A7_REG          P23
#define ZERO_REG        P24

#define NNZ_PREFETCH 128
#define ROW_PREFETCH 128

inline int __bytes_to_cachelines(int bytes)
{
    return (bytes + BSP_CONFIG_DCACHE_LINE_BYTES - 1)/BSP_CONFIG_DCACHE_LINE_BYTES;
}

void CSR_dvusmv_NxM_handler(int precision, char trans, int m, int n,
                             const double alpha,
                             const dmatCSR_t a,
                             const void * x, int x_bytes,
                             const double beta,
                             void * y, int y_bytes,
                             char enable_prefetch)
{

    if((alpha==1.0) && (beta==0.0)) {
        CSR_dvusmv_NxM_alpha1_beta0(
            precision, trans, m, n,
            a,
            x, x_bytes,
            y, y_bytes,
            enable_prefetch);
        return;
    }

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    int l_row;
    int l_col_offset;
    int l_col;
    uintptr_t l_x_ptr;
    uintptr_t l_y_ptr = (uintptr_t) y;
    uintptr_t l_val_ptr = (uintptr_t) a->val;
    
    // TODO : assert one-based indexing or change everything to 0-based

    // Load alpha and beta values in P registers
    pcvt_d_p(ALPHA_REG,	dtoraw(alpha));
    pcvt_d_p(BETA_REG, 	dtoraw(beta));

    // Iterate over line elements
    for ( l_row = 0; l_row < m; l_row++ ) {

        // acc = 0;
        pcvt_d_p(ROW_ACCU_REG, 0);

        // Iterate over column elements for the current line.
        for ( l_col_offset = ( a->ptr[l_row] -1 ) ; l_col_offset < ( a->ptr[l_row+1] -1 ) ; l_col_offset++ ) {
	
           // Get real column index
           l_col = a->ind[l_col_offset] - a->base_index;

           // Compute the address for X[l_col]
           l_x_ptr = ((uintptr_t) x) + l_col * x_bytes;
	       
           // Load X[col] in P register
           ple(X_REG, l_x_ptr, 0, EVP0);

           // Load A value in P register
           pld(A_REG, l_val_ptr, 0, EFP0);
			
           // Store the result for x[l_col] * a->val[l_col_offset] in register 
           pmul(A_REG, X_REG, A_REG, EC0);

           // acc += x[l_col] * a->val[l_col_offset] (P1)
           padd(ROW_ACCU_REG, ROW_ACCU_REG, A_REG, EC0);

           // Compute the address for A[l_col_offset]
           l_val_ptr += sizeof(double);
       }

       // acc *= alpha
       pmul(ROW_ACCU_REG, ROW_ACCU_REG, ALPHA_REG, EC0);

       // Compute the address for Y[l_row]
       l_y_ptr = ((uintptr_t) y) + l_row * y_bytes;

       // Load Y[col] in P register
       ple(Y_REG, l_y_ptr, 0, EVP1);

       // y[l_row] *= beta
       pmul(Y_REG, Y_REG, BETA_REG, EC0);

       // y[l_row] = acc * alpha + beta * y[l_row]
       padd(Y_REG, ROW_ACCU_REG, Y_REG, EC0);

       // Store Y[l_row] in memory
       pse(Y_REG, l_y_ptr, 0, EVP1);
    }
    
    VBLASPERFMONITOR_FUNCTION_END;
}

static inline __ALWAYS_INLINE__ void CSR_dvusmv_alpha1_beta0_unroll8(const uintptr_t x_ptr, int **ind_ptr, uintptr_t *val_ptr, int *elem, const int n, int x_bytes){
    int l_col_0, l_col_1, l_col_2, l_col_3,
        l_col_4, l_col_5, l_col_6, l_col_7;
    uintptr_t l_x_ptr_0, l_x_ptr_1, l_x_ptr_2, l_x_ptr_3,
        l_x_ptr_4, l_x_ptr_5, l_x_ptr_6, l_x_ptr_7;
    int _n = n;
    int *_ind_ptr = *ind_ptr;
    uintptr_t _val_ptr = *val_ptr;
    int _elem = *elem;

    pmv_p_p(ACC1_REG, ZERO_REG);
    pmv_p_p(ACC2_REG, ZERO_REG);
    pmv_p_p(ACC3_REG, ZERO_REG);
    // Optimized version
    // We try to space out instruction chains that :
    // - use the same pipeline (ADD/MUL/LSU)
    // - depends on previous result
    if (vblas_likely(_n >= 8)) {
        l_col_0  = _ind_ptr[0];
        l_col_1  = _ind_ptr[1];
        l_col_2  = _ind_ptr[2];
        l_col_3  = _ind_ptr[3];
        l_col_0 *= x_bytes;
        l_col_1 *= x_bytes;
        l_x_ptr_0  = x_ptr + l_col_0;
        l_col_2 *= x_bytes;
        ple(X0_REG, l_x_ptr_0 , 0, EVP0);
        l_x_ptr_1  = x_ptr + l_col_1;
        ple(X1_REG, l_x_ptr_1 , 0, EVP0);
        l_col_3 *= x_bytes;
        l_x_ptr_2  = x_ptr + l_col_2;
        pmv_p_p(ACC4_REG, ZERO_REG);
        ple(X2_REG, l_x_ptr_2 , 0, EVP0);
        l_x_ptr_3  = x_ptr + l_col_3;
        pmv_p_p(ACC5_REG, ZERO_REG);
        ple(X3_REG, l_x_ptr_3 , 0, EVP0);
        l_col_4  = _ind_ptr[4];
        l_col_5  = _ind_ptr[5];
        l_col_6  = _ind_ptr[6];
        l_col_7  = _ind_ptr[7];
        l_col_4 *= x_bytes;
        l_col_5 *= x_bytes;
        l_x_ptr_4  = x_ptr + l_col_4;
        l_col_6 *= x_bytes;
        ple(X4_REG, l_x_ptr_4 , 0, EVP0);
        l_x_ptr_5  = x_ptr + l_col_5;
        l_col_7 *= x_bytes;
        ple(X5_REG, l_x_ptr_5 , 0, EVP0);
        l_x_ptr_6  = x_ptr + l_col_6;
        l_x_ptr_7  = x_ptr + l_col_7;
        ple(X6_REG, l_x_ptr_6 , 0, EVP0);
        pmv_p_p(ACC6_REG, ZERO_REG);
        ple(X7_REG, l_x_ptr_7 , 0, EVP0);
        pmv_p_p(ACC7_REG, ZERO_REG);
        pld(A0_REG, _val_ptr, 0, EFP0);
        pld(A1_REG, _val_ptr, 1, EFP0);
        pld(A2_REG, _val_ptr, 2, EFP0);
        pld(A3_REG, _val_ptr, 3, EFP0);
        pld(A4_REG, _val_ptr, 4, EFP0);
        pld(A5_REG, _val_ptr, 5, EFP0);
        pld(A6_REG, _val_ptr, 6, EFP0);
        pld(A7_REG, _val_ptr, 7, EFP0);
        while(vblas_likely(_n >= 16)) {
            // We can load the next 8 elements early
            _ind_ptr += 8;
            _n-=8;
            pmul(X0_REG, X0_REG, A0_REG, EC0);
            _val_ptr += 8*sizeof(double);
            l_col_0  = _ind_ptr[0];
            // We can reuse A0 reg, we load next matrix elem
            pld(A0_REG, _val_ptr, 0, EFP0);
            pmul(X1_REG, X1_REG, A1_REG, EC0);
            l_col_1  = _ind_ptr[1];
            pld(A1_REG, _val_ptr, 1, EFP0);   
            l_col_0 *= x_bytes;            
            pmul(X2_REG, X2_REG, A2_REG, EC0);
            l_col_2  = _ind_ptr[2];       
            pld(A2_REG, _val_ptr, 2, EFP0);
            l_col_1 *= x_bytes;
            l_col_3  = _ind_ptr[3];
            pmul(X3_REG, X3_REG, A3_REG, EC0);
            pld(A3_REG, _val_ptr, 3, EFP0);
            l_col_2 *= x_bytes;
            l_col_4  = _ind_ptr[4];
            pmul(X4_REG, X4_REG, A4_REG, EC0);
            pld(A4_REG, _val_ptr, 4, EFP0);
            l_col_3 *= x_bytes;
            l_col_5  = _ind_ptr[5];
            pmul(X5_REG, X5_REG, A5_REG, EC0);
            pld(A5_REG, _val_ptr, 5, EFP0);
            l_col_4 *= x_bytes;
            l_col_6  = _ind_ptr[6];
            pmul(X6_REG, X6_REG, A6_REG, EC0);
            pld(A6_REG, _val_ptr, 6, EFP0);
            l_col_5 *= x_bytes;
            l_col_7  = _ind_ptr[7];
            pmul(X7_REG, X7_REG, A7_REG, EC0);
            pld(A7_REG, _val_ptr, 7, EFP0);
            l_x_ptr_0  = x_ptr + l_col_0;
            padd(ACC0_REG, ACC0_REG, X0_REG, EC0);
            // Here we can reuse the X0 reg, we load next
            // vector element as soon as possible
            ple(X0_REG, l_x_ptr_0 , 0, EVP0);
            l_x_ptr_1  = x_ptr + l_col_1;
            padd(ACC1_REG, ACC1_REG, X1_REG, EC0);
            ple(X1_REG, l_x_ptr_1 , 0, EVP0);
            l_x_ptr_2  = x_ptr + l_col_2;
            padd(ACC2_REG, ACC2_REG, X2_REG, EC0);
            ple(X2_REG, l_x_ptr_2 , 0, EVP0);
            l_x_ptr_3  = x_ptr + l_col_3;
            padd(ACC3_REG, ACC3_REG, X3_REG, EC0);
            ple(X3_REG, l_x_ptr_3 , 0, EVP0);
            l_x_ptr_4  = x_ptr + l_col_4;
            _elem += 8;
            padd(ACC4_REG, ACC4_REG, X4_REG, EC0);
            ple(X4_REG, l_x_ptr_4 , 0, EVP0);
            l_col_6 *= x_bytes;
            l_x_ptr_5  = x_ptr + l_col_5;
            padd(ACC5_REG, ACC5_REG, X5_REG, EC0);
            ple(X5_REG, l_x_ptr_5 , 0, EVP0);
            l_col_7 *= x_bytes;
            l_x_ptr_6  = x_ptr + l_col_6;
            padd(ACC6_REG, ACC6_REG, X6_REG, EC0);
            ple(X6_REG, l_x_ptr_6 , 0, EVP0);
            l_x_ptr_7  = x_ptr + l_col_7;
            padd(ACC7_REG, ACC7_REG, X7_REG, EC0);
            ple(X7_REG, l_x_ptr_7 , 0, EVP0);      
        } // end while(_n >= 16)
        // We have more than 8 elements left,
        // so we process 8 elements
        pmul(X0_REG, X0_REG, A0_REG, EC0);
        _val_ptr += 8*sizeof(double);
        pmul(X1_REG, X1_REG, A1_REG, EC0);
        _elem += 8;
        pmul(X2_REG, X2_REG, A2_REG, EC0); 
        _ind_ptr += 8;
        pmul(X3_REG, X3_REG, A3_REG, EC0);
        padd(ACC0_REG, ACC0_REG, X0_REG, EC0);
        pmul(X4_REG, X4_REG, A4_REG, EC0);
        padd(ACC1_REG, ACC1_REG, X1_REG, EC0);
        pmul(X5_REG, X5_REG, A5_REG, EC0);
        padd(ACC2_REG, ACC2_REG, X2_REG, EC0);
        pmul(X6_REG, X6_REG, A6_REG, EC0);
        padd(ACC3_REG, ACC3_REG, X3_REG, EC0);
        pmul(X7_REG, X7_REG, A7_REG, EC0);
        padd(ACC4_REG, ACC4_REG, X4_REG, EC0);
        padd(ACC5_REG, ACC5_REG, X5_REG, EC0);
        padd(ACC6_REG, ACC6_REG, X6_REG, EC0);
        padd(ACC7_REG, ACC7_REG, X7_REG, EC0);
        // Now we add accumulators two-by-two as we
        // won't use upper ones anymore
        _n -= 8;
        padd(ACC0_REG, ACC0_REG, ACC4_REG, EC0);
        padd(ACC1_REG, ACC1_REG, ACC5_REG, EC0);
        padd(ACC2_REG, ACC2_REG, ACC6_REG, EC0);
        padd(ACC3_REG, ACC3_REG, ACC7_REG, EC0);
    } // end if(_n >= 8)
    if(_n >= 4) {
        l_col_0  = _ind_ptr[0];
        l_col_1  = _ind_ptr[1];
        l_col_2  = _ind_ptr[2];
        l_col_3  = _ind_ptr[3];
        pld(A0_REG, _val_ptr, 0, EFP0);
        l_col_0 *= x_bytes;
        l_col_1 *= x_bytes;
        l_col_2 *= x_bytes;
        l_x_ptr_0  = x_ptr + l_col_0;
        l_col_3 *= x_bytes;
        ple(X0_REG, l_x_ptr_0 , 0, EVP0);
        l_x_ptr_1  = x_ptr + l_col_1;
        ple(X1_REG, l_x_ptr_1 , 0, EVP0);
        l_x_ptr_2  = x_ptr + l_col_2;
        pld(A1_REG, _val_ptr, 1, EFP0);
        ple(X2_REG, l_x_ptr_2 , 0, EVP0);
        l_x_ptr_3  = x_ptr + l_col_3;     
        pld(A2_REG, _val_ptr, 2, EFP0);
        ple(X3_REG, l_x_ptr_3 , 0, EVP0);
        pmul(X0_REG, X0_REG, A0_REG, EC0);
        pld(A3_REG, _val_ptr, 3, EFP0);
        _val_ptr += 4*sizeof(double);
        pmul(X1_REG, X1_REG, A1_REG, EC0);
        _ind_ptr += 4;
        _n -= 4;       
        pmul(X2_REG, X2_REG, A2_REG, EC0);
        _elem += 4;
        pmul(X3_REG, X3_REG, A3_REG, EC0);
        padd(ACC0_REG, ACC0_REG, X0_REG, EC0);
        padd(ACC1_REG, ACC1_REG, X1_REG, EC0);
        padd(ACC2_REG, ACC2_REG, X2_REG, EC0);
        padd(ACC3_REG, ACC3_REG, X3_REG, EC0);    
        // Compute the address for A[elem], ind[elem], increment elem
    } // end if(_n >= 4)
    // Add accumulators two-by-two as we won't need ACC2 and ACC3 anymore
    padd(ACC0_REG, ACC0_REG, ACC2_REG, EC0);
    padd(ACC1_REG, ACC1_REG, ACC3_REG, EC0);
    if(_n >= 2) {
        l_col_0  = _ind_ptr[0];
        l_col_1  = _ind_ptr[1];
        pld(A0_REG, _val_ptr, 0, EFP0);
        l_col_0 *= x_bytes;
        l_col_1 *= x_bytes;
        l_x_ptr_0  = x_ptr + l_col_0;
        ple(X0_REG, l_x_ptr_0 , 0, EVP0);
        l_x_ptr_1  = x_ptr + l_col_1;    
        pld(A1_REG, _val_ptr, 1, EFP0);
        ple(X1_REG, l_x_ptr_1 , 0, EVP0);
        _val_ptr += 2*sizeof(double);
        _elem += 2;
        pmul(X0_REG, X0_REG, A0_REG, EC0);
        _ind_ptr += 2;
        pmul(X1_REG, X1_REG, A1_REG, EC0);
        _n -= 2;           
        padd(ACC0_REG, ACC0_REG, X0_REG, EC0);
        padd(ACC1_REG, ACC1_REG, X1_REG, EC0);   
    } // end if(_n >= 2)
    // Add ACC0 and ACC1, ACC0 stores the result
    padd(ACC0_REG, ACC0_REG, ACC1_REG, EC0);
    if(_n >= 1) {
        l_col_0  = _ind_ptr[0];
        l_col_0 *= x_bytes;
        l_x_ptr_0  = x_ptr + l_col_0; 
        pld(A0_REG, _val_ptr, 0, EFP0);
        ple(X0_REG, l_x_ptr_0 , 0, EVP0);
        _val_ptr += 1*sizeof(double);
        _elem += 1;
        pmul(X0_REG, X0_REG, A0_REG, EC0);
        _ind_ptr += 1;
        //_n -= 1; // Not needed          
        padd(ACC0_REG, ACC0_REG, X0_REG, EC0);
    } // end if(_n >= 1) 
    // _n should be 1
    // Update pointers and elem
    *ind_ptr  = _ind_ptr;
    *val_ptr = _val_ptr;
    *elem    = _elem;
}

void CSR_dvusmv_NxM_alpha1_beta0(int precision, char trans, int m, int n,
                             const dmatCSR_t a,
                             const void * x, int x_bytes,
                             void * y, int y_bytes,
                             char enable_prefetch)
{
    int l_row;
    uintptr_t l_x_ptr;
    uintptr_t l_y_ptr = (uintptr_t) y;
    uintptr_t l_val_ptr = (uintptr_t) a->val;
    int *l_row_ptr = (int *) a->ptr;
    int *l_ind_ptr = (int *) a->ind;
    const int nnz = l_row_ptr[m]-l_row_ptr[0];

    VBLASPERFMONITOR_FUNCTION_BEGIN;

    l_x_ptr = ((uintptr_t) x) - ((a->base_index)*x_bytes); // Account for one-base indexing

#if SPVBLAS_ENABLE_HWPF
    const int PREFETCH_NNZ_PER_BLOCK = 128;
    const int PREFETCH_NNZ_BLOCKS = nnz/PREFETCH_NNZ_PER_BLOCK;
    const int PREFETCH_NNZ_REMAINING = nnz-(PREFETCH_NNZ_PER_BLOCK*PREFETCH_NNZ_BLOCKS);

    //printf("nnz : %d, nnz_per_block : %d, nnz_blocks : %d, nnz_remaining : %d\n", nnz, PREFETCH_NNZ_PER_BLOCK, PREFETCH_NNZ_BLOCKS, PREFETCH_NNZ_REMAINING);

    const int PREFETCH_ROW_PER_BLOCK = 128;

    const int HWPF_ENGINE_ROW_ID = 0;
    const int HWPF_ENGINE_IND_ID = 1;
    const int HWPF_ENGINE_VAL_ID = 2;

    hwpf_engine_params_t hwpf_param_row;
    hwpf_engine_params_t hwpf_param_ind;
    hwpf_engine_params_t hwpf_param_val;
    hwpf_engine_throttle_t hwpf_throttle_row;
    hwpf_engine_throttle_t hwpf_throttle_nnz;

    if ( enable_prefetch ) {

        hwpf_param_row.stride = 1;
        hwpf_param_row.nblocks = 1;
        hwpf_param_row.nlines = __bytes_to_cachelines(PREFETCH_ROW_PER_BLOCK*sizeof(int));

        hwpf_param_ind.stride = 1;
        hwpf_param_ind.nblocks = 1;
        hwpf_param_ind.nlines = __bytes_to_cachelines(PREFETCH_NNZ_PER_BLOCK*sizeof(int));

        hwpf_param_val.stride = 1;
        hwpf_param_val.nblocks = 1;
        hwpf_param_val.nlines = __bytes_to_cachelines(PREFETCH_NNZ_PER_BLOCK*sizeof(double));

        //printf("row nlines : %ld, ind nlines : %ld, val nlines : %ld\n", hwpf_param_row.nlines, hwpf_param_ind.nlines, hwpf_param_val.nlines);

        hwpf_throttle_row.nwait = 8; // todo: optimize
        hwpf_throttle_row.ninflight = 0;

        hwpf_throttle_nnz.nwait = 4;
        hwpf_throttle_nnz.ninflight = 0;

    //    printf("blop !\n");
        hwpf_set_params(HWPF_ENGINE_ROW_ID, &hwpf_param_row);
        hwpf_set_throttle(HWPF_ENGINE_ROW_ID, &hwpf_throttle_row);
        hwpf_set_params(HWPF_ENGINE_IND_ID, &hwpf_param_ind);
        hwpf_set_throttle(HWPF_ENGINE_IND_ID, &hwpf_throttle_nnz);
        hwpf_set_params(HWPF_ENGINE_VAL_ID, &hwpf_param_val);
        hwpf_set_throttle(HWPF_ENGINE_VAL_ID, &hwpf_throttle_nnz);
        {
            const uintptr_t _addr_row = (uintptr_t)l_row_ptr;
            hwpf_wait(HWPF_ENGINE_ROW_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_ROW_ID, _addr_row, HWPF_ENABLE);  
            const uintptr_t _addr_ind = (uintptr_t)l_ind_ptr;
            hwpf_wait(HWPF_ENGINE_IND_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_IND_ID, _addr_ind, HWPF_ENABLE);    
            const uintptr_t _addr_val = l_val_ptr;
            hwpf_wait(HWPF_ENGINE_VAL_ID);
            hwpf_set_addr_and_trigger(HWPF_ENGINE_VAL_ID, _addr_val, HWPF_ENABLE);         
        } 
        hwpf_throttle_row.nwait = 20; // todo: optimize
        hwpf_throttle_row.ninflight = 0;
        hwpf_set_throttle(HWPF_ENGINE_ROW_ID, &hwpf_throttle_row);
        hwpf_throttle_nnz.nwait = 100;
        hwpf_throttle_nnz.ninflight = 0;
        hwpf_set_throttle(HWPF_ENGINE_IND_ID, &hwpf_throttle_nnz);
        hwpf_set_throttle(HWPF_ENGINE_VAL_ID, &hwpf_throttle_nnz);
    }
#endif

    int elem = 0;
    pcvt_d_p(ZERO_REG, 0);
    int col_offset = l_row_ptr[0];
    int col_offset_next = l_row_ptr[0]; // works also for one-based indexing as it is used in difference only
    for ( l_row = 0; l_row < m; l_row++ ) {
        // init result accumulator to 0
        pmv_p_p(ACC0_REG, ZERO_REG);
#if SPVBLAS_ENABLE_HWPF
        if ( enable_prefetch ) {
            // if we are at middle of previous prefetch of row pointers
            // launch prefetch again at prefetch distance divided by 2
            if((l_row % PREFETCH_ROW_PER_BLOCK) == (PREFETCH_ROW_PER_BLOCK/2)) {
                const uintptr_t _addr_row = (uintptr_t)&l_row_ptr[PREFETCH_ROW_PER_BLOCK/2];
                hwpf_wait(HWPF_ENGINE_ROW_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_ROW_ID, _addr_row, HWPF_ENABLE);            
            }
        }
#endif
        // Compute how many columns there is in this row
        col_offset = col_offset_next;
        col_offset_next = l_row_ptr[1];
        int cols = col_offset_next-col_offset;
        // Test if the number of columns does not exceed prefetch
        //  distance divided by 2.
        // If it is the case, we may miss the prefetch trigger
        if(vblas_likely(cols < (PREFETCH_NNZ_PER_BLOCK/2))){
            // Fuzzy test to see if we will reach the middle of the current
            //  prefetch in this row
            // If yes, we need to trigger prefetch again
            if(vblas_unlikely(((elem % PREFETCH_NNZ_PER_BLOCK) <= (PREFETCH_NNZ_PER_BLOCK/2)) &&
                              (((elem+cols) % PREFETCH_NNZ_PER_BLOCK) > (PREFETCH_NNZ_PER_BLOCK/2)))) {
                // If we are at the last prefetch block, change parameters
                // to avoid going further than allocated memory
                if(elem+cols >= ((PREFETCH_NNZ_PER_BLOCK*PREFETCH_NNZ_BLOCKS)-(PREFETCH_NNZ_PER_BLOCK/2))) {
                    hwpf_param_ind.nlines = __bytes_to_cachelines(PREFETCH_NNZ_REMAINING*sizeof(int));
                    hwpf_set_params(HWPF_ENGINE_IND_ID, &hwpf_param_ind);
                    hwpf_param_val.nlines = __bytes_to_cachelines(PREFETCH_NNZ_REMAINING*sizeof(double)); 
                    hwpf_set_params(HWPF_ENGINE_VAL_ID, &hwpf_param_val);          
                }
                int dist = (((elem+PREFETCH_NNZ_PER_BLOCK-1) / PREFETCH_NNZ_PER_BLOCK) * PREFETCH_NNZ_PER_BLOCK)-elem; 
                const uintptr_t _addr_ind = (uintptr_t)&l_ind_ptr[dist];
                hwpf_wait(HWPF_ENGINE_IND_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_IND_ID, _addr_ind, HWPF_ENABLE);    
                const uintptr_t _addr_val = l_val_ptr+(sizeof(double)*dist);
                hwpf_wait(HWPF_ENGINE_VAL_ID);
                hwpf_set_addr_and_trigger(HWPF_ENGINE_VAL_ID, _addr_val, HWPF_ENABLE);    
            }
            // Execute row accumulation in unrolled fashion
            CSR_dvusmv_alpha1_beta0_unroll8(l_x_ptr, &l_ind_ptr, &l_val_ptr, &elem, cols, x_bytes);
        } else {    
            // Here we have many columns (more than prefetch distance devided
            //  by 2). So we process them block by block and check prefetch
            //  trigger after each block processing
            while(cols > 0){
                // We take the minimum between remaining columns and prefetch
                //  distance divided by 2
                int _cols = cols > PREFETCH_NNZ_PER_BLOCK/2 ? PREFETCH_NNZ_PER_BLOCK/2 : cols;
                cols -= _cols;
                // Fuzzy prefetch test, same as previous
                if(vblas_unlikely(((elem % PREFETCH_NNZ_PER_BLOCK) <= (PREFETCH_NNZ_PER_BLOCK/2)) &&
                                  (((elem+_cols) % PREFETCH_NNZ_PER_BLOCK) > (PREFETCH_NNZ_PER_BLOCK/2)))) {
                    if(elem+_cols >= ((PREFETCH_NNZ_PER_BLOCK*PREFETCH_NNZ_BLOCKS)-(PREFETCH_NNZ_PER_BLOCK/2))) {
                        hwpf_param_ind.nlines = __bytes_to_cachelines(PREFETCH_NNZ_REMAINING*sizeof(int));
                        hwpf_set_params(HWPF_ENGINE_IND_ID, &hwpf_param_ind);
                        hwpf_param_val.nlines = __bytes_to_cachelines(PREFETCH_NNZ_REMAINING*sizeof(double)); 
                        hwpf_set_params(HWPF_ENGINE_VAL_ID, &hwpf_param_val);          
                    }
                    int dist = (((elem+PREFETCH_NNZ_PER_BLOCK-1) / PREFETCH_NNZ_PER_BLOCK) * PREFETCH_NNZ_PER_BLOCK)-elem; 
                    const uintptr_t _addr_ind = (uintptr_t)&l_ind_ptr[dist];
                    hwpf_wait(HWPF_ENGINE_IND_ID);
                    hwpf_set_addr_and_trigger(HWPF_ENGINE_IND_ID, _addr_ind, HWPF_ENABLE);    
                    const uintptr_t _addr_val = l_val_ptr+(sizeof(double)*dist);
                    hwpf_wait(HWPF_ENGINE_VAL_ID);
                    hwpf_set_addr_and_trigger(HWPF_ENGINE_VAL_ID, _addr_val, HWPF_ENABLE);    
                }
                // Process _cols columns, ACC0 is not overwritten between calls
                CSR_dvusmv_alpha1_beta0_unroll8(l_x_ptr, &l_ind_ptr, &l_val_ptr, &elem, _cols, x_bytes);
            }
        }
        // Store Y[l_row] in memory
        pse(ACC0_REG, l_y_ptr, 0, EVP1);

        // Increment row and y pointers
        l_y_ptr += y_bytes;
        l_row_ptr++;
    }

    VBLASPERFMONITOR_FUNCTION_END;

}