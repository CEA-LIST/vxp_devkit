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
#ifndef __VRP_INIT_H_
#define __VRP_INIT_H_

/* VRP Reset GPIO*/

#define RST_TOP  0x1 << 0
#define RST_VRP0 0x1 << 1
#define RST_VRP1 0x1 << 2
#define RST_VRP2 0x1 << 3
#define RST_VRP3 0x1 << 4
#define RST_JTAG 0x1 << 5

//#define VRP_RELEASE_CORE0   (RST_TOP | RST_VRP0)
//#define VRP_RELEASE_TOP     RST_TOP
//#define VRP_RST_CORES       VRP_RELEASE_TOP
#define VRP_RELEASE_CORE0   0x0000000000000003
#define VRP_RELEASE_CORE1   0x0000000000000005
#define VRP_RELEASE_CORE2   0x0000000000000009
#define VRP_RELEASE_CORE3   0x0000000000000011
#define VRP_RELEASE_TOP     0x0000000000000001
#define VRP_RST_CORES       0x0000000000000001
#define VRP_RST_TOP         0x0000000000000000

#define VRP_GPIO_RST_ADDR   0x0000020240000000

/* VRP address space */

#define VRP_BASE_ADDR               0x80000000
#define VRP_BOOT_CORE0_ADDR         0x80000008
#define VRP_BOOT_CORE1_ADDR         0x80000010
#define VRP_BOOT_CORE2_ADDR         0x80000018
#define VRP_BOOT_CORE3_ADDR         0x80000020
#define VRP_IO_START_ADDR           0x80000048
#define VRP_IO_END_ADDR             0x80000050
#define VRP_CACHE_TAB_BASE0_ADDR    0x80000060
#define VRP_CACHE_TAB_MASK0_ADDR    0x80000068
#define VRP_CACHE_TAB_BASE1_ADDR    0x80000070
#define VRP_CACHE_TAB_MASK1_ADDR    0x80000078
#define VRP_CACHE_TAB_BASE2_ADDR    0x80000080
#define VRP_CACHE_TAB_MASK2_ADDR    0x80000088
#define VRP_CACHE_TAB_BASE3_ADDR    0x80000090
#define VRP_CACHE_TAB_MASK3_ADDR    0x80000098
#define VRP_AXI_STREAM_ID0_ADDR     0x80000120
#define VRP_AXI_STREAM_ID1_ADDR     0x80000128
#define VRP_AXI_STREAM_ID2_ADDR     0x80000130
#define VRP_AXI_STREAM_ID3_ADDR     0x80000138

/* VRP values */
#define VRP_IO_START        0x0000020100000000
#define VRP_IO_END          0x00000201ffffffff
#define VRP_CACHE_TAB_BASE0 0x0000020100000001
#define VRP_CACHE_TAB_MASK0 0xffffffff00000000
#define VRP_CACHE_TAB_BASE1 NULL
#define VRP_CACHE_TAB_MASK1 0xffffffffffe00000
#define VRP_EN_AXI0         0x0000000000000001
#define VRP_EN_AXI1         0x0000000000000002
#define VRP_EN_AXI2         0x0000000000000004
#define VRP_EN_AXI3         0x0000000000000008
#define VRP_AXI_STREAM_ID0  0x801
#define VRP_AXI_STREAM_ID1  0x802
#define VRP_AXI_STREAM_ID2  0x803
#define VRP_AXI_STREAM_ID3  0x804

/* VRP Files paths*/
#define VRP_WORK_FOLDER             "/usr/share/vrp"
#define VRP_FW_BINARY_PATH          "/usr/share/vrp/vrp_fw.x.bin"
#define VRP_DYN_LINK_SCRIPT_PATH    "/usr/share/vrp/scripts/vrp_dynlink.sh"
#define VRP_DYN_LINK_DONE_PATH      "/usr/share/vrp/scripts/vrp_dynlink.sh"

/* VRP Solvers values*/
#define VRP_SOLVER_STATUS_SUBMITED  0x0000000000000001
#define VRP_SOLVER_STATUS_FAILED    0x0000000000000002
#define VRP_SOLVER_STATUS_COMPLETED 0x0000000000000003

/* Functions */

int vrp_reset_platform(uint64_t base_rst_reg);
int vrp_reset_core(uint64_t base_rst_reg);
int vrp_release_core(uint64_t base_vrp_reg, uint64_t base_rst_reg);
int vrp_init_core(uint64_t base_vrp_reg, uint64_t base_rst_reg, uint64_t base_axi_reg, uint64_t vrp_fw_mem_buff);

#endif /*__VRP_INIT_H_*/
