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

#ifndef __VBLASCONFIG_HPP__
#define __VBLASCONFIG_HPP__

#include <iostream>
#include "VRPSDK/vblas_mt.h"

namespace VPFloatPackage {

    namespace VBLAS {

        struct VBLASConfig {
            // Configuration parameters for vblas_mt_init
            uint64_t nb_threads;
            uint64_t nb_rows_per_thread;

            // Configuration parameters of linear prefetcher
            uint64_t enable_prefetcher;
        };

        vblas_mt_config_t * getVBLAS_MT_Config();
        vblas_mt_vgemv_config_t * getVBLAS_MT_VGEMV_Config();

        int VBLAS_Init();

        int VBLAS_Destroy();

        VBLASConfig * VBLAS_getConfig(void);

        void VBLAS_setConfig(VBLASConfig * a_config_address);
    }
}

#endif /* __VBLASCONFIG_HPP__ */