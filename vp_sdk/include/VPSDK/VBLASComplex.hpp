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

#ifndef __VBLASCOMPLEX_HPP__
#define __VBLASCOMPLEX_HPP__

#include "VPSDK/VPComplex.hpp"
#include "Matrix/matrix.h"
#include <complex>

namespace VPFloatPackage {
    namespace VBLAS {
        /*
         *  Vector scaling - x = alpha*x
         */
        void cscal ( int n, const std::complex<float>& alpha , std::complex<float> * x , int x_inc);

        void csscal ( int n, float alpha , std::complex<float> * x, int x_inc);

        void zscal ( int n, const std::complex<double>& alpha , std::complex<double> * x, int x_inc);

        void zsscal ( int n, double alpha , std::complex<double> * x, int x_inc);

        void vcscal ( int precision , int n, const VPComplex & alpha, VPComplexArray & x);

        void vcsscal ( int precision , int n, const VPFloat & alpha, VPComplexArray & x);

        /*
         *  Vector copy - y = x
         */
        void ccopy (  int n, std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc );

        void zcopy (  int n, std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc );

        void vccopy ( int n, const VPComplexArray & x, VPComplexArray & y);

        /*
         *  Vector addition with scaling (AXPY) - y = alpha*x + y
         */
        void caxpy ( int n, const std::complex<float>& alpha , std::complex<float> * const x, int x_inc, std::complex<float> * y, int y_inc);

        void zaxpy ( int n, const std::complex<double>& alpha , std::complex<double> * const x, int x_inc, std::complex<double> * y, int y_inc);

        void vcaxpy (int precision , int n, const VPFloat & alpha, const VPComplexArray & x, VPComplexArray & y);          

        /*
         *  Scalar vector-vector multiplication (dot product) - x * y
         */
        void cdotu ( int n, std::complex<float> * const x, int x_inc, std::complex<float> * const y, int y_inc, std::complex<float> * res);
        void cdotc ( int n, std::complex<float> * const x, int x_inc, std::complex<float> * const y, int y_inc, std::complex<float> * res );
        void zdotu ( int n, std::complex<double> * const x, int x_inc, std::complex<double> * const y, int y_inc, std::complex<double> * res );
        void zdotc ( int n, std::complex<double> * const x, int x_inc, std::complex<double> * const y, int y_inc, std::complex<double> * res );
        void vcdotu (int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & r);
        void vcdotc (int precision , int n, const VPComplexArray & x, const VPComplexArray & y, VPComplex & r);        

        /*
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
         */
        void cgemv (char trans , int m, int n,
                    const std::complex<float>& alpha ,
                    std::complex<float> * const a, int lda ,
                    std::complex<float> * const x, int x_inc ,
                    const std::complex<float> & beta ,
                    std::complex<float> * y, int y_inc );

        void zgemv (char trans , int m, int n,
                    const std::complex<double>& alpha ,
                    std::complex<double> * const a, int lda ,
                    std::complex<double> * const x, int x_inc ,
                    const std::complex<double>& beta ,
                    std::complex<double> * y, int y_inc );    

        void vcgemv (int precision , char trans , int m, int n,
                     const VPComplex & alpha,
                     const matrix_t a, int lda ,
                     const VPComplexArray & x,
                     const VPComplex & beta,
                     VPComplexArray & y);

        /*****************************************************************************************************************
         * Computes the Euclidean norm of a vector. sqrt(Xdot(x))
         ****************************************************************************************************************/
        float scnrm2 (int n, std::complex<float> * const x, int x_inc, char enable_prefetch=1);

        double dznrm2 (int n, std::complex<double> * const x, int x_inc, char enable_prefetch=1);

        void vcnrm2 (int precision, int n, const VPComplexArray& x, int x_inc, VPFloat & res);
    };
};

#endif /* __VBLASCOMPLEX_HPP__ */