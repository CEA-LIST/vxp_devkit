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
#ifndef __VRP_VFIO_H__
#define __VRP_VFIO_H__

#include <linux/vfio.h>
#include <linux/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/fcntl.h>
#include <sys/mman.h>
#include <sys/eventfd.h>

#include <time.h>

#define VFIO_CONTAINER  "/dev/vfio/vfio"
#define VRP_VFIO_GROUP 	"/dev/vfio/0"      // dmesg | grep group, look for vrp address
#define VRP_DEVICE		"20100000000.vrp"   // the address for the device that is with the group

#define VRP_TRANSFER_SIZE_4KPAGES	(65536)  // 4K pages * 65536 bytes = 256MB buffer
typedef __u8 uchar;
typedef __u32 uint;
typedef __u32 u32;
typedef __u64 u64;

struct vrp_irq_efd {
    uint32_t num_irqs;
    uint32_t irq_efd[];
};

int vrp_vfio_iommu_map(int* container, int* size_to_map, uint64_t *vrp_fw_mem_buff, uint64_t *vrp_data_mem_buff, uint64_t *base_vrp_reg, uint64_t *base_rst_reg, uint64_t *base_axi_reg, struct vrp_irq_efd** vrp_efd);
int vrp_vfio_iommu_unmap(int* container, int* size_to_map, struct vrp_irq_efd** vrp_efd);

#define handle_error(msg) \
  do { perror(msg); exit(EXIT_FAILURE); } while (0)

#endif /* __VRP_VFIO_H__ */
