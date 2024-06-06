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

#ifndef __VBLAS_HPP__
#define __VBLAS_HPP__

#include "VPSDK/VPFloat.hpp"
#include "Matrix/matrix.h"

namespace VPFloatPackage {
    /*************************************************************************************************************************
     * VBLAS Namesapce
     * On a VRP platform VBLAS calls are mapped to vblas functions using vpfloat numbers.
     * On a non VRP platform VBLAS calls are mapped to openblas function using double numbers.
     ************************************************************************************************************************/
    namespace VBLAS {
        /*****************************************************************************************************************
        *  Matrix-vector multiplication
        *
        *  y = (alpha * A * x) + (beta * y)
        ****************************************************************************************************************/
        void vgemvd( int precision, char trans, int m, int n,
                            double alpha,
                            const matrix_t a,
                            const VPFloatArray & x,
                            const VPFloat & beta,
                            VPFloatArray & y);

        /*****************************************************************************************************************
         *  Vector scaling - x = alpha*x
         ****************************************************************************************************************/
        void vscal( int precision, int n, const VPFloat & alpha, VPFloatArray & x);

        /*****************************************************************************************************************
         *  Vector copy - y = x
         ****************************************************************************************************************/
        void vcopy( int n, const VPFloatArray & x, VPFloatArray & y);
        
        void vcopy_d_v( int n, const double * x, VPFloatArray & y);

        void vcopy_v_d( int n, const VPFloatArray & x, double * y);

        /*****************************************************************************************************************
         *  Vector addition with scaling (AXPY) - y = alpha*x + y
         ****************************************************************************************************************/
        void saxpy(int n, float alpha, const float *x, vpfloat_off_t x_inc, float *y, vpfloat_off_t y_inc);

        void daxpy(int n, double alpha, const double *x, vpfloat_off_t x_inc, double *y, vpfloat_off_t y_inc);

        void vaxpy( int precision, int n, const VPFloat & alpha, const VPFloatArray & x, VPFloatArray & y);

        /*****************************************************************************************************************
         *  Scalar vector-vector multiplication (dot product) - x * y
         ****************************************************************************************************************/
        float sdot(int n, const float *x, vpfloat_off_t x_inc, const float *y, vpfloat_off_t y_inc);

        double ddot(int n, const double *x, vpfloat_off_t x_inc, const double *y, vpfloat_off_t y_inc);

        void vdot( int precision, int n, const VPFloatArray & x, const VPFloatArray & y, VPFloat & res);

        void vzero(int precision, int n, VPFloatArray & x);

        /*****************************************************************************************************************
         * Computes the Euclidean norm of a vector. sqrt(Xdot(x))
         ****************************************************************************************************************/
        float snrm2 (int n, const float *x, int x_inc);

        double dnrm2 (int n, const double *x, int x_inc);

        void vnrm2 (int precision, int n, const VPFloatArray& x, VPFloat & res);
    };
};

#endif /* __VBLAS_HPP__ */
