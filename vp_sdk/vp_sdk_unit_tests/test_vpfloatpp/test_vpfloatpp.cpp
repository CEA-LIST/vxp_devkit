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
#include <VPSDK/VPFloat.hpp>

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */


int test_VPFloatArray() {
    int  l_rc = EXIT_SUCCESS;

    int l_precision=53;
    int l_exponent_size=7;
    int l_stride_size=1;
    int l_bis = l_precision + l_exponent_size + 1;
    int l_array_size = 10;

    VPFloatPackage::VPFloatComputingEnvironment::set_precision(l_bis);
    VPFloatPackage::VPFloatComputingEnvironment::set_rounding_mode(VPFloatPackage::VP_RNE);
    
    VPFloatPackage::VPFloatArray l_array(l_exponent_size, l_bis, l_stride_size, l_array_size);
    VPFloatPackage::VPFloat a(l_exponent_size, l_bis, l_stride_size);

    /* Initialize array elements */
    for (int l_index = 0 ; l_index < l_array_size; l_index++ ) {
        l_array[l_index] = l_index;
    }

    /* Check array values */
    for (int l_index = 0 ; l_index < l_array_size; l_index++ ) {
        if ( double(l_array[l_index]) != double(l_index) ) {
            std::cout << "array element " << l_index << " has value " << l_array[l_index] << " instead of " << double(l_index) << std::endl;
            l_rc = EXIT_FAILURE;
        }
    }

    /* Test that array element access and store produce a copy of the array element*/
    a = l_array[l_array_size/2];

    a += 12;

    if ( double(a) != double( ( l_array_size / 2 ) + 12 ) ) {
        std::cout << "a variable has value " << double(a) << " instead of " << double( ( l_array_size / 2 ) + 12 ) << std::endl;
        l_rc = EXIT_FAILURE;
    }

    /* Check array values again */
    for (int l_index = 0 ; l_index < l_array_size; l_index++ ) {
        if ( double(l_array[l_index]) != double(l_index) ) {
            std::cout << "array element " << l_index << " has value " << l_array[l_index] << " instead of " << double(l_index) << std::endl;
            l_rc = EXIT_FAILURE;
        }
    }

    l_array[l_array_size/2] += 12 ;

    /* Check array values again */
    for (int l_index = 0 ; l_index < l_array_size; l_index++ ) {
        if ( l_index == (l_array_size/2) ) {
            if  ( double(l_array[l_index]) != double(l_index + 12 ) ) 
            {
                std::cout << "array element " << l_index << " has value " << double(l_array[l_index]) << " instead of " << double(l_index) + 12 << std::endl;
                l_rc = EXIT_FAILURE;
                continue;
            }
        } else if ( double(l_array[l_index]) != double(l_index) ) {
            std::cout << "array element " << l_index << " has value " << double(l_array[l_index]) << " instead of " << double(l_index) << std::endl;
            l_rc = EXIT_FAILURE;
            continue;
        }       
    }

    return l_rc;
}

int main()
{
    int l_rc = EXIT_SUCCESS;

    l_rc = test_VPFloatArray(); 

    if (  test_VPFloatArray() ) {
        printf("Test VPFloatArray FAILED !\n");
        l_rc = EXIT_FAILURE;
    } else {
        printf("Test VPFloatArray OK !\n");
    }

    exit(l_rc);
}