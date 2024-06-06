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
#include "VRPSDK/vblas.h"
#include "Matrix/DENSE.h"
#include <iostream>

using namespace VPFloatPackage;

/******************************************************************************
 *  Vector scaling - x = alpha*x
 *****************************************************************************/
void VBLAS::cscal ( int n, const std::complex<float>& alpha , std::complex<float> * x , int x_inc) {
    ::cscal(n, (void *)&alpha, (void *)x, x_inc);
}

void VBLAS::csscal ( int n, float alpha , std::complex<float> * x, int x_inc) {
    ::csscal(n, alpha, (void *)x, x_inc);
}

void VBLAS::zscal ( int n, const std::complex<double>& alpha , std::complex<double> * x, int x_inc) {
    ::zscal(n, (void *)&alpha, (void *)x, x_inc);
}

void VBLAS::zsscal ( int n, double alpha , std::complex<double> * x, int x_inc) {
    ::zsscal(n, alpha, (void *)x, x_inc);
}

void VBLAS::vcscal ( int precision , int n, const VPComplex & alpha, VPComplexArray & x){
    ::vcscal(precision, n, alpha.getData(), alpha.getEnvironment(), x.getData(), x.getEnvironment());
}

void VBLAS::vcsscal ( int precision , int n, const VPFloat & alpha, VPComplexArray & x){
    ::vcsscal(precision, n, alpha.getData(), alpha.getEnvironment(), x.getData(), x.getEnvironment());
}

/******************************************************************************
 *  Vector copy - y = x
 *****************************************************************************/
void VBLAS::ccopy (  int n, std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc ) {
    ::ccopy(n, (void *)x, x_inc, (void *)y, y_inc);
}

void VBLAS::zcopy (  int n, std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc ) {
    ::zcopy(n, (void *)x, x_inc, (void *)y, y_inc);
}

void VBLAS::vccopy( int n, const VPComplexArray & x, VPComplexArray & y) {
    ::vccopy(n, x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment());
}

/******************************************************************************
 *  Vector addition with scaling (AXPY) - y = alpha*x + y
 *****************************************************************************/
void VBLAS::caxpy ( int n, const std::complex<float>& alpha , std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc) {
    ::caxpy(n, (void *)&alpha, (void *)x, x_inc, (void *)y, y_inc);
}

void VBLAS::zaxpy ( int n, const std::complex<double>& alpha , std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc) {
    ::zaxpy(n, (void *)&alpha, (void *)x, x_inc, (void *)y, y_inc);
}

void VBLAS::vcaxpy( int precision , int n, const VPFloat & alpha, const VPComplexArray & x, VPComplexArray & y){
    ::vcaxpy(precision, n, alpha.getData(), alpha.getEnvironment(), x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment()); 
}  

/******************************************************************************
 *  Scalar vector-vector multiplication (dot product) - x * y
 *****************************************************************************/
void VBLAS::cdotu ( int n, 
                    std::complex<float> * const x, int x_inc, 
                    std::complex<float> * const y, int y_inc, 
                    std::complex<float> * res) {
    ::cdotu(n, (void *)x, x_inc, (void *)y, y_inc, (void *)res);
}

void VBLAS::cdotc ( int n, 
                    std::complex<float> * const x, int x_inc, 
                    std::complex<float> * const y, int y_inc, 
                    std::complex<float> * res ) {
    ::cdotc(n, (void *)x, x_inc, (void *)y, y_inc, (void *)res);
}

void VBLAS::zdotu ( int n, 
                    std::complex<double> * const x, int x_inc, 
                    std::complex<double> * const y, int y_inc, 
                    std::complex<double> * res ) {
    ::zdotu(n, (void *)x, x_inc, (void *)y, y_inc, (void *)res);
}

void VBLAS::zdotc ( int n, 
                    std::complex<double> * const x, int x_inc, 
                    std::complex<double> * const y, int y_inc, 
                    std::complex<double> * res ) {
    ::zdotc(n, (void *)x, x_inc, (void *)y, y_inc, (void *)res);
}

void  VBLAS::vcdotu(int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & res){
    ::vcdotu(precision, n, x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment(), res.getData(), res.getEnvironment());
}

void  VBLAS::vcdotc(int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & res){
    ::vcdotc(precision, n, x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment(), res.getData(), res.getEnvironment());
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
    ::cgemv(trans, m, n, (void *)&alpha, (void *)a, lda, (void *)x, x_inc, (void *)&beta, (void *)y, y_inc);
}

void VBLAS::zgemv ( char trans , int m, int n,
                    const std::complex<double>& alpha ,
                    std::complex<double> * const a, int lda ,
                    std::complex<double> * const x, int x_inc ,
                    const std::complex<double>& beta ,
                    std::complex<double> * y, int y_inc ) {
    ::zgemv(trans, m, n, (void *)&alpha, (void *)a, lda, (void *)x, x_inc, (void *)&beta, (void *)y, y_inc);
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
            ::vcgemv(   precision, trans, m, n,
                        alpha.getData(), alpha.getEnvironment(),
                        l_a.getData(), VPFLOAT_EVP_DOUBLE,
                        lda,
                        x.getData(), x.getEnvironment(),
                        beta.getData(), beta.getEnvironment(),
                        y.getData(), y.getEnvironment());
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
    VPFloatArray l_x((float *)x, n);

    VPFloat l_res(VPFLOAT_EVP_FLOAT.es, VPFLOAT_EVP_FLOAT.bis, VPFLOAT_EVP_FLOAT.stride);

    ::vdot(VPFLOAT_EVP_FLOAT.bis, n, l_x.getData(), l_x.getEnvironment(), l_x.getData(), l_x.getEnvironment(), l_res.getData(), l_res.getEnvironment(), enable_prefetch);

    ::vsqrt(VPFLOAT_EVP_FLOAT.bis, l_res.getData(), l_res.getEnvironment(), l_res.getData(), l_res.getEnvironment());

    return float(l_res);
}

double VBLAS::dznrm2 (int n, std::complex<double> * const x, int x_inc, char enable_prefetch) {
    VPFloatArray l_x((double *)x, n);

    VPFloat l_res(VPFLOAT_EVP_DOUBLE.es, VPFLOAT_EVP_DOUBLE.bis, VPFLOAT_EVP_DOUBLE.stride);

    ::vdot(VPFLOAT_EVP_DOUBLE.bis, n, l_x.getData(), l_x.getEnvironment(), l_x.getData(), l_x.getEnvironment(), l_res.getData(), l_res.getEnvironment(), enable_prefetch);

    ::vsqrt(VPFLOAT_EVP_DOUBLE.bis, l_res.getData(), l_res.getEnvironment(), l_res.getData(), l_res.getEnvironment());

    return float(l_res);
}

void VBLAS::vcnrm2 (int precision, int n, const VPComplexArray& x, int x_inc, VPFloat & res) {
    ::vcdotu(VPFLOAT_EVP_DOUBLE.bis, n, x.getData(), x.getEnvironment(), x.getData(), x.getEnvironment(), res.getData(), res.getEnvironment());

    ::vsqrt(VPFLOAT_EVP_DOUBLE.bis, res.getData(), res.getEnvironment(), res.getData(), res.getEnvironment());
}