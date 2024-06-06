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
#include "VPSDK/VBLASConfig.hpp"
#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VMath.hpp"
#include "VRPSDK/vblas.h"
#include "VRPSDK/vblas_mt.h"
#include <stdio.h>
#include <math.h>
#include <string.h>


// #define PERF_DEBUG


#ifdef PERF_DEBUG
#include "VRPSDK/perfcounters/cpu.h"
#include <time.h>
#endif // PERF_DEBUG

extern "C" {
    #include "Matrix/matrix.h"
    #include "Matrix/BCSR.h"
    #include "Matrix/CSR.h"
    #include "Matrix/DENSE.h"
    #include "VRPSDK/spvblas.h"
    #include "bsp/bsp_config.h"
}

using namespace VPFloatPackage;
bool g_vblas_mt_initialized = false;
/*****************************************************************************************************************
*  Matrix-vector multiplication
*
*  y = (alpha * A * x) + (beta * y)
****************************************************************************************************************/
void VBLAS::vgemvd( int precision, char trans, int m, int n,
                    double alpha,
                    const matrix_t a,
                    const VPFloatArray & x,
                    const VPFloat & beta,
                    VPFloatArray & y) {
#ifdef PERF_DEBUG
    clock_t t0, t1;
    uint64_t instr0, instr1;
    uint64_t dmiss0, dmiss1;
    uint64_t imiss0, imiss1;
#endif // PERF_DEBUG

    if ( a->type_value == COMPLEX_VALUE ) {
        std::cout << __func__ << " : Matrix with complex values not supported." << std::endl;
        exit(1);
    }
		
#ifdef PERF_DEBUG		
    dmiss0 = cpu_dmiss();
    imiss0 = cpu_imiss();
    instr0 = cpu_instructions();
    t0 = clock();
#endif // PERF_DEBUG

    switch(a->type_matrix) {
    case CSR:
    case BCSR:

        if ( trans != 'N' ) {
            std::cout << __FUNCTION__ << " trans=N not supported on CSR matrix in VRP implementation." << std::endl;
            return;
        }

        ::dvusmv(   precision, trans,
                    alpha,
                    a,
                    x.getData(), x.getEnvironment(),
                    (double)beta,
                    y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

        break;

    case DENSE:

    if ( BSP_CONFIG_NCPUS > 1 && getVBLAS_MT_Config()->max_threads > 0 ) {
        ::vgemv_mt(VPFloatPackage::VBLAS::getVBLAS_MT_VGEMV_Config(), 
            precision, trans, m, n,
            (void *)&alpha, VPFLOAT_EVP_DOUBLE,
            ((dmatDENSE_t)(a->matrix->repr))->val, VPFLOAT_EVP_DOUBLE,
            a->lda,
            x.getData(), x.getEnvironment(),
            beta.getData(), beta.getEnvironment(),
            y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);
    } else {
        ::vgemv(precision, trans, m, n,
            (void *)&alpha, VPFLOAT_EVP_DOUBLE,
            ((dmatDENSE_t)(a->matrix->repr))->val, VPFLOAT_EVP_DOUBLE,
            a->lda,
            x.getData(), x.getEnvironment(),
            beta.getData(), beta.getEnvironment(),
            y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);
    }

        break;
    default:
        printf("%s : Matrix type %ld is not supported\n", __func__, a->matrix->type_id);
        std::cout << "Matrix type " << a->matrix->type_id << " not supported." << std::endl;
        break;
    };

#ifdef PERF_DEBUG
    t1 = clock();
    dmiss1 = cpu_dmiss();
    imiss1 = cpu_imiss();
    instr1 = cpu_instructions();

    double ipc = ((double)(instr1-instr0)/(t1-t0));

    std::cout << __func__ << " instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << std::endl;
#endif // PERF_DEBUG

}

/*****************************************************************************************************************
 *  Vector scaling - x = alpha*x
 ****************************************************************************************************************/
void VBLAS::vscal( int precision, int n, const VPFloat & alpha, VPFloatArray & x) {
#ifdef PERF_DEBUG
    clock_t t0, t1;
    uint64_t instr0, instr1;
    uint64_t dmiss0, dmiss1;
    uint64_t imiss0, imiss1;
    dmiss0 = cpu_dmiss();
    imiss0 = cpu_imiss();
    instr0 = cpu_instructions();
    t0 = clock();
#endif // PERF_DEBUG

    ::vscal(precision, n, alpha.getData(), alpha.getEnvironment(), x.getData(), x.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

#ifdef PERF_DEBUG
    t1 = clock();
    dmiss1 = cpu_dmiss();
    imiss1 = cpu_imiss();
    instr1 = cpu_instructions();

    double ipc = ((double)(instr1-instr0)/(t1-t0));

    std::cout << __func__ << " instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << std::endl;
#endif // PERF_DEBUG 
}

/*****************************************************************************************************************
 *  Vector copy - y = x
 ****************************************************************************************************************/
void VBLAS::vcopy( int n, const VPFloatArray & x, VPFloatArray & y) {
#ifdef PERF_DEBUG
    clock_t t0, t1;
    uint64_t instr0, instr1;
    uint64_t dmiss0, dmiss1;
    uint64_t imiss0, imiss1;
    dmiss0 = cpu_dmiss();
    imiss0 = cpu_imiss();
    instr0 = cpu_instructions();
    t0 = clock();
#endif // PERF_DEBUG

    ::vcopy(n, x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

#ifdef PERF_DEBUG
    t1 = clock();
    dmiss1 = cpu_dmiss();
    imiss1 = cpu_imiss();
    instr1 = cpu_instructions();

    double ipc = ((double)(instr1-instr0)/(t1-t0));

    std::cout << __func__ << " instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << std::endl;
#endif // PERF_DEBUG 
}

void VBLAS::vcopy_d_v( int n, const double * x, VPFloatArray & y) {

    ::vcopy_d_v(n, x, 1, y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

}

void VBLAS::vcopy_v_d( int n, const VPFloatArray & x, double * y) {

    ::vcopy_v_d(n, x.getData(), x.getEnvironment(), y, 1, VBLAS_getConfig()->enable_prefetcher);

}


/*****************************************************************************************************************
 *  Vector addition with scaling (AXPY) - y = alpha*x + y
 ****************************************************************************************************************/
void VBLAS::saxpy(int n, float alpha, const float *x, vpfloat_off_t x_inc, float *y, vpfloat_off_t y_inc) {

    ::vaxpy( VPFLOAT_EVP_FLOAT.bis, n,
                (void *)&alpha, VPFLOAT_EVP_FLOAT,
                x, VPFLOAT_EVP_FLOAT,
                y, VPFLOAT_EVP_FLOAT,
                VBLAS_getConfig()->enable_prefetcher);

}

void VBLAS::daxpy(int n, double alpha, const double *x, vpfloat_off_t x_inc, double *y, vpfloat_off_t y_inc) { 

    ::vaxpy( VPFLOAT_EVP_DOUBLE.bis, n,
                (void *)&alpha, VPFLOAT_EVP_DOUBLE,
                x, VPFLOAT_EVP_DOUBLE,
                y, VPFLOAT_EVP_DOUBLE,
                VBLAS_getConfig()->enable_prefetcher);

}

void VBLAS::vaxpy( int precision, int n, const VPFloat & alpha, const VPFloatArray & x, VPFloatArray & y) {
#ifdef PERF_DEBUG
    clock_t t0, t1;
    uint64_t instr0, instr1;
    uint64_t dmiss0, dmiss1;
    uint64_t imiss0, imiss1;
    dmiss0 = cpu_dmiss();
    imiss0 = cpu_imiss();
    instr0 = cpu_instructions();
    t0 = clock();
#endif // PERF_DEBUG  

::vaxpy(precision, n, alpha.getData(), alpha.getEnvironment(), x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

#ifdef PERF_DEBUG    
    t1 = clock();
    dmiss1 = cpu_dmiss();
    imiss1 = cpu_imiss();
    instr1 = cpu_instructions();

    double ipc = ((double)(instr1-instr0)/(t1-t0));

    std::cout << __func__ << " instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << std::endl;
#endif // PERF_DEBUG
}

/*****************************************************************************************************************
 *  Scalar vector-vector multiplication (dot product) - x * y
 ****************************************************************************************************************/
float VBLAS::sdot(int n, const float *x, vpfloat_off_t x_inc, const float *y, vpfloat_off_t y_inc) {               
    VPFloatPackage::VPFloat l_result((float)0.0);

    ::vdot( VPFLOAT_EVP_FLOAT.bis, 
            n,
            x, VPFLOAT_EVP_FLOAT,
            y, VPFLOAT_EVP_FLOAT,
            l_result.getData(), l_result.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

    return float(l_result);
}

double VBLAS::ddot(int n, const double *x, vpfloat_off_t x_inc, const double *y, vpfloat_off_t y_inc) {
    VPFloatPackage::VPFloat l_result((double)0.0);

    ::vdot( VPFLOAT_EVP_DOUBLE.bis, 
            n,
            x, VPFLOAT_EVP_DOUBLE,
            y, VPFLOAT_EVP_DOUBLE,
            l_result.getData(), l_result.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

    return double(l_result);
}

void VBLAS::vdot( int precision, int n, const VPFloatArray & x, const VPFloatArray & y, VPFloat & res) {
#ifdef PERF_DEBUG
    clock_t t0, t1;
    uint64_t instr0, instr1;
    uint64_t dmiss0, dmiss1;
    uint64_t imiss0, imiss1;
    dmiss0 = cpu_dmiss();
    imiss0 = cpu_imiss();
    instr0 = cpu_instructions();
    t0 = clock();
#endif // PERF_DEBUG 
    
    ::vdot(precision, n, x.getData(), x.getEnvironment(), y.getData(), y.getEnvironment(), res.getData(), res.getEnvironment(), VBLAS_getConfig()->enable_prefetcher);

#ifdef PERF_DEBUG
    t1 = clock();
    dmiss1 = cpu_dmiss();
    imiss1 = cpu_imiss();
    instr1 = cpu_instructions();

    double ipc = ((double)(instr1-instr0)/(t1-t0));

    std::cout << __func__ << " instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << std::endl;
#endif // PERF_DEBUG    
}

void VBLAS::vzero(int precision, int n, VPFloatArray & x) {

    memset(x.getData(), 0, n*VPFLOAT_SIZEOF(x.getEnvironment()));

}

/*****************************************************************************************************************
 * Computes the Euclidean norm of a vector. sqrt(Xdot(x))
 ****************************************************************************************************************/
float VBLAS::snrm2 (int n, const float *x, int x_inc) {

    return sqrtf(::sdot(n, x, x_inc, x, x_inc, VBLAS_getConfig()->enable_prefetcher));

}

double VBLAS::dnrm2 (int n, const double *x, int x_inc) {

    return sqrt(::ddot(n, x, x_inc, x, x_inc, VBLAS_getConfig()->enable_prefetcher));

}

void VBLAS::vnrm2 (int precision, int n, const VPFloatArray& x, VPFloat & res) {

    vdot(precision, n, x, x, res);

    vsqrt(precision, res.getData(), res.getEnvironment(), res.getData(), res.getEnvironment());

}
