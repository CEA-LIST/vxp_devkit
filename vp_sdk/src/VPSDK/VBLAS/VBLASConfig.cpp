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
#include <string.h>
using namespace VPFloatPackage::VBLAS;

static VBLASConfig g_vblas_config = {   .nb_threads = 1,
                                        .nb_rows_per_thread = 64,
                                        .enable_prefetcher = 1};

vblas_mt_config_t * VPFloatPackage::VBLAS::getVBLAS_MT_Config() {
    static vblas_mt_config_t vblas_mt_config;

    vblas_mt_config.max_threads = ::g_vblas_config.nb_threads;

    return &vblas_mt_config;
}

vblas_mt_vgemv_config_t * VPFloatPackage::VBLAS::getVBLAS_MT_VGEMV_Config() {
    static vblas_mt_vgemv_config_t vblas_mt_vgemv_config;

    vblas_mt_vgemv_config.min_rows_per_job = ::g_vblas_config.nb_rows_per_thread;
    vblas_mt_vgemv_config.vblas_mt = getVBLAS_MT_Config();

    return &vblas_mt_vgemv_config;
}

VBLASConfig* VPFloatPackage::VBLAS::VBLAS_getConfig(void) {
    return &::g_vblas_config;
}

void VPFloatPackage::VBLAS::VBLAS_setConfig(VBLASConfig* a_config_address) {
    memcpy(&g_vblas_config, a_config_address, sizeof(VBLASConfig));

    std::cout << "VBLASConfig => nb_thread : " << g_vblas_config.nb_threads << std::endl;
    std::cout << "VBLASConfig => nb_rows_per_thread : " << g_vblas_config.nb_rows_per_thread << std::endl;
    std::cout << "VBLASConfig => enable_prefetcher : " << g_vblas_config.enable_prefetcher << std::endl;
}