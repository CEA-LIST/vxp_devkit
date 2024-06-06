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
 * Authors       : Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#include "VPSDK/VBLAS.hpp"
#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VMath.hpp"
#include <mpfr.h>
#include <math.h>
#include <stdio.h>
#include "Matrix/matrix.h"
#include "Matrix/BCSR.h"
#include "Matrix/CSR.h"
#include "Matrix/DENSE.h"
using namespace VPFloatPackage;

/*****************************************************************************************************************
*  Matrix-vector multiplication
*
*  y = (alpha * A * x) + (beta * y)
****************************************************************************************************************/
void vgemvdBCSR(int precision, char trans, int m, int n,
                double alpha,
                const dmatBCSR_t a_bcsr,
                int lda,
                const VPFloatArray & x,
                const VPFloat & beta,
                VPFloatArray & y, 
                VPFloatArray & acc,
                int start_row_number) {
    int l_nb_elements_per_block = a_bcsr->row_block_size * a_bcsr->col_block_size;
    int l_block_index = 0;
    int l_block_val_index = 0;

    // Processing of row of blocks
    for ( int l_row_block = 0 ; l_row_block < a_bcsr->num_block_rows; l_row_block++ ) {
        
        int l_real_start_row_index = l_row_block * a_bcsr->row_block_size + start_row_number;

        // Processing a block
        for (int l_block = a_bcsr->bptr[l_row_block]; l_block < a_bcsr->bptr[l_row_block+1]; l_block++) {

            int l_real_start_col_index = a_bcsr->bind[l_block_index];

            // Process each block line
            for ( int l_row_in_block = 0; l_row_in_block < a_bcsr->row_block_size; l_row_in_block++ ) {

                int l_real_row_index = l_real_start_row_index + l_row_in_block;

                // Process each column in a block line
                for ( int l_col_in_block = 0; l_col_in_block < a_bcsr->col_block_size; l_col_in_block++ ) {

                    int l_element_in_block_offset = l_row_in_block * a_bcsr->col_block_size + l_col_in_block;

                    int l_real_col_index = l_real_start_col_index + l_col_in_block;
                    
                    acc[l_real_row_index] += x[l_real_col_index] * a_bcsr->bval[l_block_val_index + l_element_in_block_offset];
                }

            }

            l_block_val_index += l_nb_elements_per_block;
            l_block_index++;

        }
    }

    if ( a_bcsr->num_rows_leftover > 0 ) {
        vgemvdBCSR( precision, trans, m, n,
                    alpha,
                    a_bcsr->leftover,
                    lda, 
                    x, 
                    beta, 
                    y,
                    acc,
                    start_row_number + ( a_bcsr->num_block_rows * a_bcsr->row_block_size));
    }
}

void VBLAS::vgemvd( int precision, char trans, int m, int n,
                    double alpha,
                    const matrix_t a,
                    const VPFloatArray & x,
                    const VPFloat & beta,
                    VPFloatArray & y) {

    if ( a->type_value == COMPLEX_VALUE ) {
        std::cout << __func__ << " : Matrix with complex values not supported." << std::endl;
    }

    switch(a->type_matrix) {

        case CSR: {
            if ( trans != 'N' ) {
                std::cout << __FUNCTION__ << " trans=N not supported on CSR matrix in MPFR implementation." << std::endl;
                return;
            }

            VPFloat res(VPFloatComputingEnvironment::get_temporary_var_environment().es,
                        VPFloatComputingEnvironment::get_temporary_var_environment().bis,
                        VPFloatComputingEnvironment::get_temporary_var_environment().stride);
            int j, i,k;

            dmatCSR_t l_csr = (dmatCSR_t)a->matrix->repr;

            // i est l'indice de ligne
            for (i=0; i<m; i++) {
                res = 0.0;

                // k est l'indice de colonne
                for (k=(l_csr->ptr[i])-1; k<(l_csr->ptr[i+1])-1; k++) {
                    // k est l'indice de ligne
                    j = (l_csr->ind[k])-1; // tjrs le decalage magique

                    res +=  x[j] * l_csr->val[k];
                }

                y[i] = res * alpha + beta * y[i];
            }
        }; break;

        case BCSR: {
            dmatBCSR_t l_bcsr = (dmatBCSR_t)a->matrix->repr;
            vpfloat_evp_t y_env = y.getEnvironment();

            if ( trans != 'N' ) {
                std::cout << __FUNCTION__ << " trans=N not supported on BCSR matrix in MPFR implementation." << std::endl;
                return;
            }
            

            VPFloatArray l_acc(y_env.es, y_env.bis, y_env.stride, y.nbElements());

            for ( int i = 0 ; i < y.nbElements(); i++ ) {
                l_acc[i] = double(0.0);
            }

            vgemvdBCSR( precision, trans, m, n,
                        alpha,
                        l_bcsr,
                        a->lda,
                        x,
                        beta,
                        y, 
                        l_acc, 0);

            for (int l_row_index = 0 ; l_row_index < m; l_row_index++){
                y[l_row_index] =  y[l_row_index] * beta + l_acc[l_row_index] * alpha;
            }
        }; break;
            
        case DENSE: {
            int i,k;
            double * l_dense = (double *)(((dmatDENSE_t)(a->matrix->repr))->val);
            VPFloat res(VPFloatComputingEnvironment::get_temporary_var_environment().es,
                        VPFloatComputingEnvironment::get_temporary_var_environment().bis,
                        VPFloatComputingEnvironment::get_temporary_var_environment().stride);

            for (i=0; i<m; i++) {
                res=0.0;
                for (k=0; k<n; k++) {
                    if ( trans == 'N' ) {
                        res += x[k] * l_dense[(i*a->lda)+k];
                    } else {
                        res += x[k] * l_dense[(k * a->lda) + i];
                    }
                }
                y[i] = res * alpha + beta * y[i];
            }

        }; break;
        default: 
            std::cout << "Matrix type " << a->matrix->type_id << " not supported." << std::endl;
            break;
    };
}

/*****************************************************************************************************************
 *  Vector scaling - x = alpha*x
 ****************************************************************************************************************/
void VBLAS::vscal( int precision, int n, const VPFloat & alpha, VPFloatArray & x) {
    int i;

    for (i=0; i<n; i++) {
        x[i]=alpha * x[i]; 
    }
}

/*****************************************************************************************************************
 *  Vector copy - y = x
 ****************************************************************************************************************/
void VBLAS::vcopy( int n, const VPFloatArray & x, VPFloatArray & y) {
    int i;

    for (i=0; i<n; i++) {
        y[i]=x[i];
    }
}

void VBLAS::vcopy_d_v( int n, const double * x, VPFloatArray & y) {
    int i;

    for (i=0; i<n; i++) {
        y[i]=x[i];
    }
}

void VBLAS::vcopy_v_d( int n, const VPFloatArray & x, double * y) {
    int i;

    for (i=0; i<n; i++) {
        y[i]=x[i];
    }
}


/*****************************************************************************************************************
 *  Vector addition with scaling (AXPY) - y = alpha*x + y
 ****************************************************************************************************************/
void VBLAS::saxpy(int n, float alpha, const float *x, vpfloat_off_t x_inc, float *y, vpfloat_off_t y_inc) {
    int i;

    for (i=0; i<n; i++) {
        y[i*y_inc]=alpha * x[i*x_inc]+y[i*y_inc]; 
    }
}

void VBLAS::daxpy(int n, double alpha, const double *x, vpfloat_off_t x_inc, double *y, vpfloat_off_t y_inc) {
    int i;

    for (i=0; i<n; i++) {
        y[i*y_inc]=alpha * x[i*x_inc]+y[i*y_inc]; 
    }
}

void VBLAS::vaxpy( int precision, int n, const VPFloat & alpha, const VPFloatArray & x, VPFloatArray & y) {
    int i;

    for (i=0; i<n; i++) {
        y[i] = alpha * x[i]+y[i]; 
    }    

}

/*****************************************************************************************************************
 *  Scalar vector-vector multiplication (dot product) - x * y
 ****************************************************************************************************************/
float VBLAS::sdot(int n, const float *x, vpfloat_off_t x_inc, const float *y, vpfloat_off_t y_inc) {               
    float ftemp=0.0;
    int i;

    for (i=0; i<n; i++) {
        ftemp += y[i*y_inc]*x[i*x_inc];
    }
    
    return(ftemp);
}

double VBLAS::ddot(int n, const double *x, vpfloat_off_t x_inc, const double *y, vpfloat_off_t y_inc) {
    double dtemp=0.0;
    int i;

    for (i=0; i<n; i++) {
        dtemp += y[i*y_inc]*x[i*x_inc];
    }

    return(dtemp);
}

void VBLAS::vdot( int precision, int n, const VPFloatArray & x, const VPFloatArray & y, VPFloat & res) {
    int i;
    res = 0.0;
    for (i=0; i<n; i++) {
        res += y[i]*x[i];
    }
}

void VBLAS::vzero(int precision, int n, VPFloatArray & x) {
    int l_index;

    for (l_index = 0 ;  l_index < x.nbElements(); l_index++) {
        x[l_index] = 0.0;
    }
}

/*****************************************************************************************************************
 * Computes the Euclidean norm of a vector. sqrt(Xdot(x))
 ****************************************************************************************************************/
float VBLAS::snrm2 (int n, const float *x, int x_inc) {
    return sqrtf(VBLAS::sdot(n, x, x_inc, x, x_inc));
}

double VBLAS::dnrm2 (int n, const double *x, int x_inc) {
    return sqrt(VBLAS::ddot(n, x, x_inc, x, x_inc));
}

void VBLAS::vnrm2 (int precision, int n, const VPFloatArray& x, VPFloat & res) {
    vdot(precision, n, x, x, res);
    const VPFloat & l_const_res = res;
    res = VMath::vsqrt(l_const_res);
}