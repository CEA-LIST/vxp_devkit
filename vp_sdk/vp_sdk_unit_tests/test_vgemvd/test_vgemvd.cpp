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

#include <stdio.h>

#include <iostream>
#include <iomanip>

#include "Matrix/matrix.h"
#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VBLAS.hpp"
#include "OSKIHelper.hpp"

void vgemvd(int a_precision, const matrix_t a_matrix, const VPFloatPackage::VPFloatArray & a_x, double a_alpha, VPFloatPackage::VPFloat & a_beta, int a_n, double * a_y_val) {
    VPFloatPackage::VPFloatArray l_y(a_y_val, a_n);  
 
    VPFloatPackage::VBLAS::vgemvd(a_precision, 'N', a_n, a_n, a_alpha, a_matrix, a_x, a_beta, l_y);

    for (int i = 0 ; i < a_n; i++ ){
        a_y_val[i] = double(l_y[i]);
    }
}

void BCSRvgemvd(int a_precision, const matrix_t a_matrix, const VPFloatPackage::VPFloatArray & a_x, double a_alpha, VPFloatPackage::VPFloat & a_beta, int a_n, double * a_y_val) {
    VPFloatPackage::VPFloatArray l_y(a_y_val, a_n);  

    oski_matrix_wrapper_t l_oski_csr_matrix = VPFloatPackage::OSKIHelper::fromCSRMatrix(a_matrix);

    matrix_t l_bcsr_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_csr_matrix, 3, 3);

    if (l_bcsr_matrix->matrix->repr == NULL ) {
        printf("%s FAILED!\n", __FUNCTION__);
        return;
    }

    vgemvd(a_precision, l_bcsr_matrix, a_x, a_alpha, a_beta, a_n, a_y_val);

}

void DENSEvgemvd(int a_precision, const matrix_t a_matrix, const VPFloatPackage::VPFloatArray & a_x, double a_alpha, VPFloatPackage::VPFloat & a_beta, int a_n, double * a_y_val) {
    VPFloatPackage::VPFloatArray l_y(a_y_val, a_n);  

    oski_matrix_wrapper_t l_oski_csr_matrix = VPFloatPackage::OSKIHelper::fromCSRMatrix(a_matrix);

    matrix_t l_dense_matrix = VPFloatPackage::OSKIHelper::toDense(l_oski_csr_matrix, false, 0);

    if (l_dense_matrix->matrix->repr == NULL ) {
        printf("%s FAILED!\n", __FUNCTION__);
        return;
    }

    vgemvd(a_precision, l_dense_matrix, a_x, a_alpha, a_beta, a_n, a_y_val);
}

int main(int argc, char *argv[])
{
    int l_n = 10;
    int l_precision = 512;
    short l_exponent_size = 7;
    short l_stride_size = 1;
    short l_bis=l_precision+l_exponent_size+l_stride_size;

    int l_row_ptr[]={1,3,6,7,9,10,10,12,12,13,17};
    int l_col_ind[]={1,10,2,3,7,5,4,10,5,4,7,2,1,6,9,10};
    double l_val[]={10,220,20,40,90,60,70,210,80,120,110,100,200,160,180,190};

    double l_alpha = 10.0;
    VPFloatPackage::VPFloat l_beta(2.0);
    

    double * l_x_val = (double *)malloc(sizeof(double) * l_n);
    for (int i = 0 ; i < l_n; i++ ){
        l_x_val[i] = double(i+1);
    }

    VPFloatPackage::VPFloatArray l_x(l_x_val, l_n);  
    matrix_t l_matrix = VPFloatPackage::OSKIHelper::toMatrix(VPFloatPackage::OSKIHelper::buildCSR(l_n, l_n, l_row_ptr, l_col_ind, l_val, 1), false);
    // matrix_t l_matrix = buildCSR(10, 10, l_row_ptr, l_col_ind, l_val, 1);

    VPFloatPackage::VPFloatComputingEnvironment::set_precision(l_precision);
    VPFloatPackage::VPFloatComputingEnvironment::set_tempory_var_environment(l_exponent_size, l_bis, l_stride_size);

    double * l_y_csr_val = (double *)malloc(sizeof(double) * l_n);
    double * l_y_bcsr_val = (double *)malloc(sizeof(double) * l_n);
    double * l_y_dense_val = (double *)malloc(sizeof(double) * l_n);

    for (int i = 0 ; i < l_n; i++ ){
        l_y_csr_val[i] = l_y_bcsr_val[i] = l_y_dense_val[i] = double(l_n - (i));
    }

    vgemvd(l_precision, l_matrix, l_x, l_alpha, l_beta, l_n, l_y_csr_val);

    DENSEvgemvd(l_precision, l_matrix, l_x, l_alpha, l_beta, l_n, l_y_dense_val);

    BCSRvgemvd(l_precision, l_matrix, l_x, l_alpha, l_beta, l_n, l_y_bcsr_val);

    bool l_diff_detected = false;

    for (int i = 0 ; i < l_n; i++ ){
        // std::cout  << "i : " << i << " - dense : " << std::setw(6) << l_y_dense_val[i] << " - csr : " << std::setw(6) << l_y_csr_val[i] << " - bcsr : " << std::setw(6) << l_y_bcsr_val[i] << std::endl;

        if ( l_y_dense_val[i] != l_y_csr_val[i] || l_y_dense_val[i] != l_y_bcsr_val[i] ) {
            l_diff_detected = true;
        }
    }

    if ( l_diff_detected ) {
        std::cout << "ERROR : Difference detected!" << std::endl;
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
        exit(0);
    }
}


