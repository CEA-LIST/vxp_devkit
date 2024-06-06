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
 * Authors       : Riccardo ALIDORI
 * Creation Date : August, 2022
 * Description   : 
 **/
/******************************************************************************
 *
 *  @file   vrp_init.c : VRP handling routines
 *  @author Riccardo ALIDORI
 *
******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "vrp_init.h"

int vrp_reset_platform(uint64_t base_rst_reg) {
    uint64_t addr;
    uint64_t data;

    // Reset the platform
    addr = base_rst_reg + 0x00;
    data = VRP_RST_TOP;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
    return 0;
}

int vrp_reset_core(uint64_t base_rst_reg) {
    uint64_t addr;
    uint64_t data;

    // Reset the VRP cores
    addr = base_rst_reg + 0x00;
    data = VRP_RST_CORES;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
    return 0;
}

int vrp_release_core(uint64_t base_vrp_reg, uint64_t base_rst_reg) {
    uint64_t addr;
    uint64_t data;

    // Release VRP TOP and VRP 0
    addr = base_rst_reg + 0x00;
    data = VRP_RELEASE_CORE0 | VRP_RELEASE_CORE1;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Enable AXI
    addr = base_vrp_reg + VRP_BASE_ADDR;
    data = VRP_EN_AXI0 | VRP_EN_AXI1;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
    return 0;
}

int vrp_init_core(uint64_t base_vrp_reg, uint64_t base_rst_reg, uint64_t base_axi_reg, uint64_t vrp_fw_mem_buff) {
    uint64_t addr;
    uint64_t data;

    uint64_t boot_addr = (uint64_t)vrp_fw_mem_buff + 0x200000;

    // Force AXI CACHE to 4b'1xxx and AXI PROT to 3b'x1x
    addr = base_axi_reg + 0x00;
    data = 0x28;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Release Top Reset
    addr = base_rst_reg + 0x00;
    data = VRP_RELEASE_TOP;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

	// Setting boot address CORE 0
    addr = base_vrp_reg + VRP_BOOT_CORE0_ADDR;
    data = boot_addr;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
    //
	// Setting boot address CORE 1
    addr = base_vrp_reg + VRP_BOOT_CORE1_ADDR;
    data = boot_addr;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

	// Setting boot address CORE 2
    addr = base_vrp_reg + VRP_BOOT_CORE2_ADDR;
    data = boot_addr;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

	// Setting boot address CORE 3
    addr = base_vrp_reg + VRP_BOOT_CORE3_ADDR;
    data = boot_addr;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting IO Segment base
    addr = base_vrp_reg + VRP_IO_START_ADDR;
    data = VRP_IO_START;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting IO Segment end
    addr = base_vrp_reg + VRP_IO_END_ADDR;
    data = VRP_IO_END;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting Cacheable table 0 base + EN
    addr = base_vrp_reg + VRP_CACHE_TAB_BASE0_ADDR;
    data = VRP_CACHE_TAB_BASE0;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting Cacheable table 0 mask
    addr = base_vrp_reg + VRP_CACHE_TAB_MASK0_ADDR;
    data = VRP_CACHE_TAB_MASK0;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting Cacheable table 1 base + EN
    addr = base_vrp_reg + VRP_CACHE_TAB_BASE1_ADDR;
    data = vrp_fw_mem_buff+0x1;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
    //
    // Setting Cacheable table 0 mask
    addr = base_vrp_reg + VRP_CACHE_TAB_MASK1_ADDR;
    data = VRP_CACHE_TAB_MASK1;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting AXI Stream ID0 for the IOMMU
    addr = base_vrp_reg + VRP_AXI_STREAM_ID0_ADDR;
    data = VRP_AXI_STREAM_ID0;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting AXI Stream ID1 for the IOMMU
    addr = base_vrp_reg + VRP_AXI_STREAM_ID1_ADDR;
    data = VRP_AXI_STREAM_ID1;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting AXI Stream ID2 for the IOMMU
    addr = base_vrp_reg + VRP_AXI_STREAM_ID2_ADDR;
    data = VRP_AXI_STREAM_ID2;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // Setting AXI Stream ID3 for the IOMMU
    addr = base_vrp_reg + VRP_AXI_STREAM_ID3_ADDR;
    data = VRP_AXI_STREAM_ID3;
    printf("%s : 0x%lx = 0x%lx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    return 0;
}
