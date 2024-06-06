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
#include <stdlib.h>

#include <iostream>
#include <iomanip>

#include "OSKIHelper.hpp"
#include "Matrix/matrix.h"

int test_get(matrix_t a_matrix, double * a_val) {
    double l_value = get(a_matrix, -1, 6);
    if (! isnan(l_value)) {
        std::cout << "FAIL" << " - 1) get with invalid dimension does not return NAN " << l_value<< std::endl;    
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    l_value = get(a_matrix, 6, -10);
    if (! isnan(l_value)) {
        std::cout << "FAIL" << " - 2) get with invalid dimension does not return NAN" << std::endl;    
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    l_value = get(a_matrix, -10, -10);
    if (! isnan(l_value)) {
        std::cout << "FAIL" << " - 3) get with invalid dimension does not return NAN" << std::endl;    
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    l_value = get(a_matrix, 7, 7);
    if ( l_value != a_val[9] ) {
        std::cout << "FAIL" << " - 4) get with valid dimension does not return right value " << a_val[9] << " " << l_value << std::endl;    
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
    }
    
    l_value = get(a_matrix, 10, 6);
    if ( l_value != a_val[12] ) {
        std::cout << "FAIL" << " - 5) get with valid dimension does not return right value " << a_val[12] << " " << l_value << std::endl;    
        exit(1);
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    return 0;
}

void test_for_diag_matrix_build(int a_base_index) {
    int n=10;
    double *myval;

    /* creer la matrice */
    myval = (double *)malloc(n*sizeof(double));

    for (int l_value_index=0; l_value_index < n; l_value_index++) {
        myval[l_value_index]=double(l_value_index+1); 
    }

    matrix_t l_iM = buildDiagCSR(n, myval, a_base_index);

    std::cout << "=== Diagonal CSR manualy created matrix with base_index " << a_base_index << " display test ====" << std::endl;
    displayMatrix(l_iM, 3);

    std::cout << "=== Diagonal CSR manualy created matrix with base_index " << a_base_index << " get test ====" << std::endl;
    double l_i_previous_value = get(l_iM, a_base_index, a_base_index);

    for (int l_value_index = 1 + a_base_index; l_value_index < n + a_base_index; l_value_index++) {
        double l_i_value = get(l_iM, l_value_index, l_value_index);

        if ( ( l_i_value - l_i_previous_value ) != 1 ) {
            std::cout << "FAIL" << " Diag values (" << l_value_index << ", " << l_value_index << ") not equal to " << l_value_index << " " << l_i_value << " " << l_i_previous_value << std::endl;
            exit(1);
        }

        l_i_previous_value = l_i_value;
    }

    std::cout << "SUCCESS" << std::endl;
}

int main(int argc, char *argv[])
{
    int l_rc = 0;
    int l_row_ptr[]={1,3,6,6,8,9,9,11,11,12,16};
    int l_col_ind[]={1,10,2,3,7,4,10,5,4,7,2,1,6,9,10};
    double l_val[]={10,220,20,40,90,70,210,80,120,110,100,200,160,180,190};
    int l_base_index = 1;

    const char* l_matrix_repo_path = getenv("MATRIX_REPO_PATH");

    matrix_t l_matrix = buildCSR(10, 10, l_row_ptr, l_col_ind, l_val, l_base_index);

    oski_matrix_wrapper_t l_simple_oski_CSR_matrix = VPFloatPackage::OSKIHelper::fromCSRMatrix(l_matrix);

    matrix_t l_dense_matrix = VPFloatPackage::OSKIHelper::toDense(l_simple_oski_CSR_matrix);

    matrix_t l_bcsr_matrix_extracted_from_oski = VPFloatPackage::OSKIHelper::toBCSR(l_simple_oski_CSR_matrix, 4 , 4);
    
    std::cout << "=== Manualy built CSR Matrix display test ====" << std::endl;
    displayMatrix(l_matrix, 3);

    std::cout << "=== DENSE Matrix display test ====" << std::endl;
    displayMatrix(l_dense_matrix, 3);

    std::cout << "=== BCSR Matrix display test ====" << std::endl;
    displayMatrix(l_bcsr_matrix_extracted_from_oski, 3);

    std::cout << "=== DENSE Matrix get test ====" << std::endl;
    l_rc = test_get(l_dense_matrix, l_val);

    std::cout << "=== CSR Matrix get test ====" << std::endl;
    l_rc = test_get(l_matrix, l_val);

    std::cout << "=== BCSR Matrix get test ====" << std::endl;
    l_rc = test_get(l_bcsr_matrix_extracted_from_oski, l_val);

    std::cout << "=== BIG CSR matrix display test ====" << std::endl;
    char * l_matrix_file_path = (char *)malloc(sizeof(char) * ( strlen(l_matrix_repo_path) + strlen("/bcsstk01.mtx") + 1 ));
    sprintf(l_matrix_file_path, "%s/bcsstk01.mtx", l_matrix_repo_path);

    oski_matrix_wrapper_t l_oski_csr_matrix = VPFloatPackage::OSKIHelper::loadFromFile(l_matrix_file_path);

    matrix_t l_csr_matrix_extracted_from_oski = VPFloatPackage::OSKIHelper::toMatrix(l_oski_csr_matrix);
    displayMatrix(l_csr_matrix_extracted_from_oski, 1);

    test_for_diag_matrix_build(l_base_index);
    test_for_diag_matrix_build(0);

    std::cout << "SUCCESS" << std::endl;
    exit(l_rc);
}