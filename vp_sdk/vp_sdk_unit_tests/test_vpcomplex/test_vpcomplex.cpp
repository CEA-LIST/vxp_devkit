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

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */


using namespace VPFloatPackage;

int test_add(uint16_t es, uint16_t ms, uint16_t stride, VPComplex& a, VPComplex& b) {
    VPComplex l_complex_add(es,ms,stride);

    l_complex_add = a + b;

    if ( double(l_complex_add.real()) != double(a.real() + b.real()) || double(l_complex_add.imag()) != double(a.imag() + b.imag()) ) {
        std::cout << "l_complex_add : real : " << double(l_complex_add.real()) << " - imag : " << double(l_complex_add.imag()) << std::endl;
        std::cout << "expected      : real : " << double(a.real() + b.real()) << " - imag : " << double(a.imag() + b.imag()) << std::endl;        
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int test_sub(uint16_t es, uint16_t ms, uint16_t stride, VPComplex& a, VPComplex& b) {
    VPComplex l_complex_sub(es,ms,stride);

    l_complex_sub = a - b;

    if ( double(l_complex_sub.real()) != double(a.real() - b.real()) || double(l_complex_sub.imag()) != double(a.imag() - b.imag()) ) {
        std::cout << "l_complex_sub : real : " << double(l_complex_sub.real()) << " - imag : " << double(l_complex_sub.imag()) << std::endl;
        std::cout << "expected      : real : " << double(a.real() - b.real()) << " - imag : " << double(a.imag() - b.imag()) << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


int test_mul(uint16_t es, uint16_t ms, uint16_t stride, VPComplex& a, VPComplex& b) {
    VPComplex l_complex_mul(es,ms,stride);

    l_complex_mul = a * b;

    if ( double(l_complex_mul.real()) != double(a.real() * b.real() - a.imag() * b.imag()) || double(l_complex_mul.imag()) != double(a.real() * b.imag() + a.imag() * b.real()) ) {
        std::cout << "l_complex_mul : real : " << double(l_complex_mul.real()) << " - imag : " << double(l_complex_mul.imag()) << std::endl;
        std::cout << "expected      : real : " << double(a.real() * b.real() - a.imag() * b.imag()) << " - imag : " << double(a.real() * b.imag() + a.imag() * b.real()) << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int test_conjugate(uint16_t es, uint16_t ms, uint16_t stride, VPComplex& a) {
    VPComplex l_a_conjugate(es,ms,stride);

    l_a_conjugate = a.conjugate();

    if ( double(l_a_conjugate.real()) != double(a.real()) || double(l_a_conjugate.imag()) != -double(a.imag()) ) {
        std::cout << "l_a_conjugate : real : " << double(l_a_conjugate.real()) << " - imag : " << double(l_a_conjugate.imag()) << std::endl;
        std::cout << "expected      : real : " << double(a.real()) << " - imag : " << -double(a.imag()) << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int test_div(uint16_t es, uint16_t ms, uint16_t stride, VPComplex& a, VPComplex& b) {
    VPComplex l_complex_div(es,ms,stride);
    VPComplex l_numerator(es,ms,stride);
    VPComplex l_denominator(es,ms,stride);
    VPFloat l_div_result_real(es,ms,stride);
    VPFloat l_div_result_imag(es,ms,stride);
    VPComplex l_b_conjugate(es,ms,stride);

    VPComplex l_local_div_result(es,ms,stride);

    l_b_conjugate = b.conjugate();

    l_complex_div = a / b;

    l_numerator = a * l_b_conjugate;
    l_denominator = b * l_b_conjugate;

    l_local_div_result.real(l_numerator.real() / l_denominator.real());
    l_local_div_result.imag(l_numerator.imag() / l_denominator.real());

    if ( l_complex_div != l_local_div_result ) {
        std::cout << "l_complex_div : real : " << double(l_complex_div.real()) << " - imag : " << double(l_complex_div.imag()) << std::endl;
        std::cout << "expected      : real : " << double(l_div_result_real) << " - imag : " << double(l_div_result_imag) << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

int test_VPComplex() {

    int l_rc = EXIT_SUCCESS;

    uint16_t es = 11;
    uint16_t ms = 53;
    uint16_t stride = 1;

    VPFloatComputingEnvironment::set_precision(53);
    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_tempory_var_environment(es, ms, 1);

    VPComplex l_a(es,ms,stride);
    VPComplex l_b(es,ms,stride);

    l_a.real(3.0);
    l_a.imag(-6.0);

    l_b.real(1.0);
    l_b.imag(-5.0);

    /*
     * Check Add operation
     */    
    if ( test_add(es, ms, stride, l_a, l_b) ){
        printf("Test complex addition FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test complex addition OK !\n");
    }
    /*
     * Check Conjugate computation
     */
    if ( test_conjugate(es, ms, stride, l_a) ) {
        printf("Test complex conjugate FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test complex conjugate OK !\n");
    }

    /*
     * Check substraction operation
     */
    if ( test_sub(es, ms, stride, l_a, l_b) ) {
        printf("Test complex substraction FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test complex substraction OK !\n");
    }

    /*
     * Check multiplication operation
     */
    if ( test_mul(es, ms, stride, l_a, l_b) ) {
        printf("Test complex multiplication FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test complex multiplication OK !\n");
    }

    /*
     * Check Division operation
     */
    if ( test_div(es, ms, stride, l_a, l_b) ) {
        printf("Test complex division FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test complex division OK !\n");
    }

    return l_rc;
}

int main()
{
    int l_rc = EXIT_SUCCESS;

    l_rc = test_VPComplex(); 

    exit(l_rc);
}