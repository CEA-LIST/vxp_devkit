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

#include "VPSDK/VBLASComplex.hpp"
#include "VPSDK/VMath.hpp"
#include "Matrix/DENSE.h"
#include <iostream>

using namespace VPFloatPackage;

/******************************************************************************
 *  Vector scaling - x = alpha*x
 *****************************************************************************/
void VBLAS::cscal ( int n, const std::complex<float>& alpha , std::complex<float> * x , int x_inc) {
    for (int x_index=0; x_index < n; x_index += x_inc) {
        x[x_index] = x[x_index] * alpha;
    }
}

void VBLAS::csscal ( int n, float alpha , std::complex<float> * x, int x_inc) {
    for (int x_index=0; x_index < n; x_index += x_inc) {
        x[x_index] = x[x_index] * alpha;
    }
}

void VBLAS::zscal ( int n, const std::complex<double>& alpha , std::complex<double> * x, int x_inc) {
    for (int x_index=0; x_index < n; x_index += x_inc) {
        x[x_index] = x[x_index] * alpha;
    }
}

void VBLAS::zsscal ( int n, double alpha , std::complex<double> * x, int x_inc) {
    for (int x_index=0; x_index < n; x_index += x_inc) {
        x[x_index] = x[x_index] * alpha;
    }
}

void VBLAS::vcscal ( int precision , int n, const VPComplex & alpha, VPComplexArray & x){
    for (int i=0; i<n; i++) {
        x[i]=x[i] * alpha; 
    }
}

void VBLAS::vcsscal ( int precision , int n, const VPFloat & alpha, VPComplexArray & x){
    for (int i=0; i<n; i++) {
        x[i]= x[i] * alpha; 
    }
}

/******************************************************************************
 *  Vector copy - y = x
 *****************************************************************************/
void VBLAS::ccopy (  int n, std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc ) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            y[y_index] = x[x_index];
            y_index += y_inc;
    }
}

void VBLAS::zcopy (  int n, std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc ) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            y[y_index] = x[x_index];
            y_index += y_inc;
    }
}

void VBLAS::vccopy( int n, const VPComplexArray & x, VPComplexArray & y) {
    for (int i=0; i<n; i++) {
        y[i]=x[i]; 
    }
}

/******************************************************************************
 *  Vector addition with scaling (AXPY) - y = alpha*x + y
 *****************************************************************************/
void VBLAS::caxpy ( int n, const std::complex<float>& alpha , std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            y[y_index] += alpha * x[x_index];
            y_index += y_inc;
    }
}

void VBLAS::zaxpy ( int n, const std::complex<double>& alpha , std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            y[y_index] += alpha * x[x_index];
            y_index += y_inc;
    }
}

void VBLAS::vcaxpy( int precision , int n, const VPFloat & alpha, const VPComplexArray & x, VPComplexArray & y){
    for (int i=0; i<n; i++) {
        y[i] += x[i] * alpha; 
    }
}  

/******************************************************************************
 *  Scalar vector-vector multiplication (dot product) - x * y
 *****************************************************************************/
void VBLAS::cdotu ( int n, 
                    std::complex<float> * const x, int x_inc, 
                    std::complex<float> * const y, int y_inc, 
                    std::complex<float> * res) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            res[x_index] = x[x_index] * y[y_index];
            y_index += y_inc;
    }           
}

void VBLAS::cdotc ( int n, 
                    std::complex<float> * const x, int x_inc, 
                    std::complex<float> * const y, int y_inc, 
                    std::complex<float> * res ) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            res[x_index] = x[x_index] * y[y_index];
            y_index += y_inc;
    }  
}

void VBLAS::zdotu ( int n, 
                    std::complex<double> * const x, int x_inc, 
                    std::complex<double> * const y, int y_inc, 
                    std::complex<double> * res ) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            res[x_index] = x[x_index] * y[y_index];
            y_index += y_inc;
    } 
}

void VBLAS::zdotc ( int n, 
                    std::complex<double> * const x, int x_inc, 
                    std::complex<double> * const y, int y_inc, 
                    std::complex<double> * res ) {
    int y_index = 0;
    for ( int x_index = 0 ; x_index < n  && y_index < n; x_index += x_inc ) {
            res[x_index] = x[x_index] * y[y_index];
            y_index += y_inc;
    } 
}

void VBLAS::vcdotu(int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & r){
    int i;
    
    r.real(0.0);
    r.imag(0.0);

    for (i=0; i<n; i++) {
        r += x[i]*y[i];
    }
}

void  VBLAS::vcdotc(int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & r){
    int i;
    
    r.real(0.0);
    r.imag(0.0);

    for (i=0; i<n; i++) {
        r += x[i].conjugate()*y[i];
    }
}  

/******************************************************************************
 *  Matrix-vector multiplication
 *
 *  y = (alpha * A * x) + (beta * y)
 *
 *  @param [in]      precision     Working precision of computations
 *  @param [in]      trans         'N' no transposition of A. 'T' transpose A
 *  @param [in]      m             Number of rows in A
 *  @param [in]      n             Number of columns in A
 *  @param [in]      alpha         Pointer to scalar
 *  @param [in]      alpha_evp     Memory environment of alpha scalar
 *  @param [in]      a             Pointer to the A matrix
 *  @param [in]      a_evp         Memory environment of A matrix
 *  @param [in]      lda           Leading dimension of A (allows spacing between rows)
 *  @param [in]      x             Pointer to vector X
 *  @param [in]      x_evp         Memory environment of X vector
 *  @param [in]      b             Pointer to scalar
 *  @param [in]      b_evp         Memory environment of beta scalar
 *  @param [in,out]  y             Pointer to vector Y
 *  @param [in]      y_evp         Memory environment of Y vector
 *****************************************************************************/
void VBLAS::cgemv ( char trans , int m, int n,
                    const std::complex<float>& alpha ,
                    std::complex<float> * const a, int lda ,
                    std::complex<float> * const x, int x_inc ,
                    const std::complex<float>& beta ,
                    std::complex<float> * y, int y_inc ) {
    int i,k;
    
    std::complex<float> res;;

    for (i=0; i<m; i++) {
        res.real(0.0);
        res.imag(0.0);
        for (k=0; k<n; k++) {
            if ( trans != 'N' ) {
                res += x[k] * a[(k*m)+i];
            } else {
                res += x[k] * a[(i*n)+k];
            }
        }
        y[i] = res * alpha + beta*y[i];
    }  
}

void VBLAS::zgemv ( char trans , int m, int n,
                    const std::complex<double>& alpha ,
                    std::complex<double> * const a, int lda ,
                    std::complex<double> * const x, int x_inc ,
                    const std::complex<double>& beta ,
                    std::complex<double> * y, int y_inc ) {
    int i,k;
    
    std::complex<double> res;;

    for (i=0; i<m; i++) {
        res.real(0.0);
        res.imag(0.0);
        for (k=0; k<n; k++) {
            if ( trans != 'N' ) {
                res += x[k] * a[(k*m)+i];
            } else {
                res += x[k] * a[(i*n)+k];
            }
        }
        y[i] = res * alpha + beta*y[i];
    }      
}

void vcgemv_DENSE (int precision , char trans , int m, int n,
                    const VPComplex & alpha,
                    const VPComplexArray & a, int lda ,
                    const VPComplexArray & x,
                    const VPComplex & beta,
                    VPComplexArray & y) {
    int i,k;
    
    VPComplex res(  VPFloatComputingEnvironment::get_temporary_var_environment().es,
                    VPFloatComputingEnvironment::get_temporary_var_environment().bis,
                    VPFloatComputingEnvironment::get_temporary_var_environment().stride);

    for (i=0; i<m; i++) {
        res.real(0.0);
        res.imag(0.0);
        for (k=0; k<n; k++) {
            if ( trans != 'N' ) {
                res += x[k] * a[(k*m)+i];
            } else {
                res += x[k] * a[(i*n)+k];
            }
        }
        y[i] = res * alpha + beta*y[i];
    }    
}

void VBLAS::vcgemv (int precision , char trans , int m, int n,
                    const VPComplex & alpha,
                    const matrix_t a, int lda ,
                    const VPComplexArray & x,
                    const VPComplex & beta,
                    VPComplexArray & y) {
    switch(a->type_matrix) {  
        case DENSE: {
            VPComplexArray l_a((double *)(((zmatDENSE_t)(a->matrix->repr))->val), m * n);
            vcgemv_DENSE(precision, trans, m, n, alpha, l_a, lda, x, beta, y);
        }; break;
        default:
            std::cout << "Matrix type " << a->matrix->type_id << " not supported." << std::endl;
            break;
    }
}   


/*****************************************************************************************************************
 * Computes the Euclidean norm of a vector. sqrt(Xdot(x))
 ****************************************************************************************************************/
float VBLAS::scnrm2 (int n, std::complex<float> * const x, int x_inc, char enable_prefetch) {
    double l_result = 0.0;
    for ( int i = 0; i < n; i += x_inc) {
        l_result += pow(std::norm(x[i]), 2.0);
    }

    return std::sqrt(l_result);
}

double VBLAS::dznrm2 (int n, std::complex<double> * const x, int x_inc, char enable_prefetch) {
    double l_result = 0.0;
    for ( int i = 0; i < n; i += x_inc) {
        l_result += pow(std::norm(x[i]), 2.0);
    }

    return std::sqrt(l_result);
}

void VBLAS::vcnrm2 (int precision, int n, const VPComplexArray & x, int x_inc, VPFloat & res) {
    VPFloat l_result(res.getEnvironment().es, res.getEnvironment().bis, res.getEnvironment().stride);

    l_result = 0.0;

    for ( int i = 0; i < n; i += x_inc) {
        VPFloat l_x_norm(x[i].getEnvironment().es, x[i].getEnvironment().bis, x[i].getEnvironment().stride);
        l_x_norm = x[i].norm();
        l_result += l_x_norm * l_x_norm;
    }

    res = VMath::vsqrt(l_result);
}

