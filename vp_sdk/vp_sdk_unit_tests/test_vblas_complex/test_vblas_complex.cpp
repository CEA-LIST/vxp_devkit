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

#include <iostream>
#include <cstdlib>
#include <VPSDK/VPComplex.hpp>
#include <VPSDK/VBLASComplex.hpp>
#include <VPSDK/VMath.hpp>
#include "OSKIHelper.hpp"
#include "Matrix/DENSE.h"

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */


using namespace VPFloatPackage;

VPComplex myadd(const VPComplex & a, const VPComplex & b) {
    //(ac-bd)+(ad+bc)i;
    VPComplex l_result(a.getEnvironment().es, a.getEnvironment().bis, a.getEnvironment().stride);

    l_result.real(a.real() + b.real());
    l_result.imag(a.imag() + b.imag());

    return l_result;
}

VPComplex mymult(const VPComplex & a, const VPFloat & b) {
    //(ac-bd)+(ad+bc)i;
    VPComplex l_result(a.getEnvironment().es, a.getEnvironment().bis, a.getEnvironment().stride);

    l_result.real(a.real()*b);
    l_result.imag(a.imag()*b);

    return l_result;
}

VPComplex mymult(const VPComplex & a, const VPComplex & b) {
    //(ac-bd)+(ad+bc)i;
    VPComplex l_result(a.getEnvironment().es, a.getEnvironment().bis, a.getEnvironment().stride);

    l_result.real(a.real()*b.real() - a.imag()*b.imag());
    l_result.imag(a.real()*b.imag() + a.imag()*b.real());

    return l_result;
}

int test_scnrm2() {
    int l_rc = EXIT_SUCCESS;

    int n = 10;
    std::complex<float> l_x[n];
    float l_norm2 = 0.0;
    float l_check_norm = 0.0;

    for ( int i = 0 ; i < n ; i++ ) {
        l_x[i].real(( i + 1 ) * 3);
        l_x[i].imag(( i + 3 ) * 3);

        l_check_norm += std::pow(std::norm(l_x[i]), 2.0);
    }

    l_check_norm = sqrtf(l_check_norm);

    l_norm2 = VBLAS::scnrm2 (n, l_x, 1);

    if ( l_check_norm != l_norm2 ) {
        std::cout << "computed value: " << l_norm2 << " -- expected value:" << l_check_norm << std::endl;
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_dznrm2() {
    int l_rc = EXIT_SUCCESS;

    int n = 10;
    std::complex<double> l_x[n];
    double l_norm2 = 0.0;
    double l_check_norm = 0.0;

    for ( int i = 0 ; i < n ; i++ ) {
        l_x[i].real(( i + 1 ) * 3);
        l_x[i].imag(( i + 3 ) * 3);

        l_check_norm += std::pow(std::norm(l_x[i]), 2.0);
    }

    l_check_norm = sqrt(l_check_norm);

    l_norm2 = VBLAS::dznrm2 (n, l_x, 1);

    if ( l_check_norm != l_norm2 ) {
        std::cout << "computed value: " << l_norm2 << " -- expected value:" << l_check_norm << std::endl;
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_vcnrm2() {
    int l_rc = EXIT_SUCCESS;

    int n = 10;
    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    VPComplexArray l_x(es, ms, stride, n);
    VPFloat l_norm2(es, ms, stride);
    VPFloat l_check_norm(es, ms, stride);

    for ( int i = 0 ; i < n ; i++ ) {
        l_x[i].real(( i + 1 ) * 3);
        l_x[i].imag(( i + 3 ) * 3);

        VPFloat l_x_norm(es, ms, stride);
        l_x_norm = l_x[i].norm();

        l_check_norm += l_x_norm * l_x_norm;
    }

    const VPFloat & l_const_check_norm = l_check_norm;
    l_check_norm = VMath::vsqrt(l_const_check_norm);

    VBLAS::vcnrm2 (ms, n, l_x, 1, l_norm2);

    if ( l_check_norm != l_norm2 ) {
        std::cout << "computed value: " << l_norm2 << " -- expected value:" << l_check_norm << std::endl;
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_VBLAS_vcgemv() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    int n = 5;
    int m = n;
    int l_diff_count = 0;
    char l_trans = 'N';
    const int l_mtx_filename_size=30;
    char l_mtx_filename[l_mtx_filename_size];

    VPFloatComputingEnvironment::set_precision(ms);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);
    
    VPComplex l_alpha(es, ms, stride);
    VPComplex l_beta(es, ms, stride);
    VPComplexArray l_x_vector(es, ms, stride, n);
    VPComplexArray l_y_vector(es, ms, stride, n);
    VPComplexArray l_backup_vector(es, ms, stride, n);
    VPComplex l_tmp_result(es, ms, stride);
    VPComplex l_row_result(es, ms, stride);

    // Initialize alpha and beta
    l_alpha.real(6.0);
    l_alpha.imag(8.0);

    l_beta.real(18.0);
    l_beta.imag(16.0);

    // Initialize X and Y vector
    for ( int i = 0 ; i< n; i++ ){
        l_x_vector[i].real(double(i));
        l_x_vector[i].imag(double(-i));
        l_y_vector[i].real(double(i) * 11);
        l_y_vector[i].imag(double(-i) * 11);
    }

    snprintf(l_mtx_filename, l_mtx_filename_size, "./test_vblas_complex/test.mtx");
    oski_matrix_wrapper_t l_oski_sparse_input_matrix = VPFloatPackage::OSKIHelper::loadFromFile(l_mtx_filename);
    matrix_t l_a = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix);

    VBLAS::vccopy(n, l_y_vector, l_backup_vector);
    VBLAS::vcgemv(ms, l_trans, m, n, l_alpha, l_a, 1, l_x_vector, l_beta, l_y_vector);

    VPComplexArray l_a_data((double *)(((zmatDENSE_t)(l_a->matrix->repr))->val), m * n);

    for (int i=0; i<m; i++) {
        l_tmp_result.real(0.0);
        l_tmp_result.imag(0.0);
        l_row_result.real(0.0);
        l_row_result.imag(0.0);

        for (int k=0; k<n; k++) {
            if ( l_trans != 'N' ) {
                l_tmp_result = mymult(l_x_vector[k], l_a_data[(k*m)+i]);
            } else {
                l_tmp_result = mymult(l_x_vector[k], l_a_data[(i*n)+k]);
            }
            l_row_result = myadd(l_row_result, l_tmp_result);
        }
        l_tmp_result = mymult(l_beta, l_backup_vector[i]);
        l_row_result = mymult(l_row_result, l_alpha);
        l_row_result = myadd(l_row_result, l_tmp_result);

        if ( l_row_result != l_y_vector[i] ) {
            std::cout << l_y_vector[i] << " -- " << l_row_result << std::endl;
            l_diff_count++;
        }
    }  

    if ( l_diff_count ) {
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_VBLAS_vcdotX() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    int n = 10;

    VPFloatComputingEnvironment::set_precision(ms);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);

    VPComplexArray l_x_vector(es, ms, stride, n);
    VPComplexArray l_y_vector(es, ms, stride, n);
    VPComplex l_dotresult(es, ms, stride);
    VPComplex l_result(es, ms, stride);
    VPComplex l_x_conjugate(es, ms, stride);

    l_dotresult.real(0.0);
    l_dotresult.imag(0.0);
    l_result.real(0.0);
    l_result.imag(0.0);

    // Initialize X vector
    for ( int i = 0 ; i< n; i++ ){
        l_x_vector[i].real(double(i));
        l_x_vector[i].imag(double(-i));
        l_y_vector[i].real(double(i) * 11);
        l_y_vector[i].imag(double(-i) * 11);
    }
    
    VBLAS::vcdotu(ms , n, l_x_vector, l_y_vector, l_dotresult);

    for ( int i = 0 ; i < n; i++ ){
        VPComplex l_tmp = mymult(l_x_vector[i], l_y_vector[i]);
        l_result = myadd(l_result, l_tmp);
    }

    if ( l_dotresult != l_result ) {
        l_rc = EXIT_FAILURE;
    }

    VBLAS::vcdotc(ms , n, l_x_vector, l_y_vector, l_dotresult);

    l_result.real(0.0);
    l_result.imag(0.0);

    for ( int i = 0 ; i < n; i++ ){
        l_x_conjugate.real(l_x_vector[i].real());
        l_x_conjugate.imag(-l_x_vector[i].imag());
        VPComplex l_tmp = mymult(l_x_conjugate, l_y_vector[i]);
        l_result = myadd(l_result, l_tmp);
    }

    if ( l_dotresult != l_result ) {
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_VBLAS_vcaxpy() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    int n = 10;
    int l_diff_count = 0;

    VPFloatComputingEnvironment::set_precision(ms);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);

    VPFloat l_alpha(es, ms, stride);
    VPComplexArray l_x_vector(es, ms, stride, n);
    VPComplexArray l_y_vector(es, ms, stride, n);
    VPComplexArray l_backup_vector(es, ms, stride, n);

    l_alpha = 11.0;

    // Initialize X vector
    for ( int i = 0 ; i< n; i++ ){
        l_x_vector[i].real(double(i));
        l_x_vector[i].imag(double(-i));
        l_y_vector[i].real(double(i) * 11);
        l_y_vector[i].imag(double(-i) * 11);
    }
    
    VBLAS::vccopy(n, l_y_vector, l_backup_vector);
    VBLAS::vcaxpy(ms , n, l_alpha, l_x_vector, l_y_vector);

    for ( int i = 0 ; i < n; i++ ){
        VPComplex l_tmp = mymult(l_x_vector[i], l_alpha);
        VPComplex l_result = myadd(l_tmp, l_backup_vector[i]);

        if ( l_y_vector[i] != l_result ) {
            std::cout << l_y_vector[i] << " -- " << l_result << std::endl;
            l_diff_count++;
        }
    }

    if ( l_diff_count ) {
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_VBLAS_vcsscal() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    int n = 10;
    int l_diff_count = 0;

    VPFloatComputingEnvironment::set_precision(ms);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);

    VPFloat l_alpha(es, ms, stride);
    VPComplexArray l_x_vector(es, ms, stride, n);
    VPComplexArray l_y_vector(es, ms, stride, n);

    l_alpha = 11.0;

    // Initialize X vector
    for ( int i = 0 ; i< n; i++ ){
        l_x_vector[i].real(double(i));
        l_x_vector[i].imag(double(-i));
    }
    
    VBLAS::vccopy(n, l_x_vector, l_y_vector);
    VBLAS::vcsscal(ms , n, l_alpha, l_x_vector);

    for ( int i = 0 ; i < n; i++ ){
        VPComplex l_result=mymult(l_y_vector[i], l_alpha);
        
        if ( l_x_vector[i] != l_result ) {
            std::cout << l_x_vector[i] << " -- " << l_result << std::endl;
            l_diff_count++;
        }
    }

    if ( l_diff_count ) {
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int test_VBLAS_vcscal() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;
    int n = 10;
    int l_diff_count = 0;

    VPFloatComputingEnvironment::set_precision(ms);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);

    VPComplex l_alpha(es, ms, stride);
    VPComplexArray l_x_vector(es, ms, stride, n);
    VPComplexArray l_y_vector(es, ms, stride, n);

    l_alpha.real(6.0);
    l_alpha.imag(8.0);

    // Initialize X vector
    for ( int i = 0 ; i< n; i++ ){
        l_x_vector[i].real(double(i));
        l_x_vector[i].imag(double(-i));
    }
    
    VBLAS::vccopy(n, l_x_vector, l_y_vector);
    VBLAS::vcscal(ms , n, l_alpha, l_x_vector);

    for ( int i = 0 ; i < n; i++ ){
        VPComplex l_result=mymult(l_y_vector[i], l_alpha);
        
        if ( l_x_vector[i] != l_result ) {
            std::cout << l_x_vector[i] << " -- " << l_result << std::endl;
            l_diff_count++;
        }
    }

    if ( l_diff_count ) {
        l_rc = EXIT_FAILURE;
    }

    return l_rc;
}

int main()
{
    int l_rc = EXIT_SUCCESS;

    if ( test_VBLAS_vcscal() ) {
        printf("Test VBLAS complex vcscal FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcscal OK !\n");
    }

    if ( test_VBLAS_vcsscal() ) {
        printf("Test VBLAS complex vcsscal FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcsscal OK !\n");
    }

    if ( test_VBLAS_vcaxpy() ) {
        printf("Test VBLAS complex vcaxpy FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcaxpy OK !\n");
    }

    if ( test_VBLAS_vcdotX() ) {
        printf("Test VBLAS complex vcdotX FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcdotX OK !\n");
    }

    if ( test_VBLAS_vcgemv() ) {
        printf("Test VBLAS complex vcgemv FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcgemv OK !\n");
    }

    if ( test_scnrm2() ) {
        printf("Test VBLAS complex scnrm2 FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex scnrm2 OK !\n");
    }

    if ( test_dznrm2() ) {
        printf("Test VBLAS complex dznrm2 FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex dznrm2 OK !\n");
    }

    if ( test_vcnrm2() ) {
        printf("Test VBLAS complex vcnrm2 FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VBLAS complex vcnrm2 OK !\n");
    }

    exit(l_rc);
}