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
 * Creation Date : October, 2023
 * Description   :
 **/


#include "VPSDK/VBLASConfig.hpp"
#include "VRPSDK/vblas_mt.h"
#include "common/cache.h"

using namespace VPFloatPackage::VBLAS;

int VPFloatPackage::VBLAS::VBLAS_Init() {
    int l_rc = EXIT_SUCCESS;
    vblas_mt_config_t * l_vblas_mt_config = getVBLAS_MT_Config();

    if ( l_vblas_mt_config->max_threads > 0 ) {
        cpu_dcache_invalidate();

        if ( ::vblas_mt_init(l_vblas_mt_config) < 0 ) {
            std::cout << __FUNCTION__ << " Fail initializing vblas_mt environment!" << std::endl;
            l_rc = EXIT_FAILURE;
        } else {
            std::cout << __FUNCTION__ << " vblas_mt environment initialized!" << std::endl;        
        }
    } else {
        std::cout << __FUNCTION__ << " No vblas_environment destruction needed!" << std::endl;
    }

    return l_rc;
}

int VPFloatPackage::VBLAS::VBLAS_Destroy() {
    int l_rc = EXIT_SUCCESS;
    vblas_mt_config_t * l_vblas_mt_config = getVBLAS_MT_Config();

    if ( l_vblas_mt_config->max_threads > 0 ) {

        if ( ::vblas_mt_destroy(getVBLAS_MT_Config()) < 0 ) {
            std::cout << __FUNCTION__ << " Fail destroying vblas_mt environment!" << std::endl;
            l_rc = EXIT_FAILURE;
        } else {
            std::cout << __FUNCTION__ << " vblas_mt environment destroyed!" << std::endl;
        }

    } else {
        std::cout << __FUNCTION__ << " No vblas_environment destruction needed!" << std::endl;
    }

    return l_rc;
}
