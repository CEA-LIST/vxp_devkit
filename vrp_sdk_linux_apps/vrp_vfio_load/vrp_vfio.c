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
 *  @file   vrp_vfio.c : Maps and Unmaps the VRP in the iommu (SMMU)
 *  @author Riccardo ALIDORI
 *
******************************************************************************/

#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include "vrp_vfio.h"

int vrp_vfio_iommu_map(int* container, int* size_to_map, uint64_t *vrp_fw_mem_buff, uint64_t *vrp_data_mem_buff, uint64_t *base_vrp_reg, uint64_t *base_rst_reg, uint64_t *base_axi_reg, struct vrp_irq_efd** vrp_efd) {
    int l_rc = 0;

    int group, device;
	unsigned int i;

	struct vfio_group_status group_status = { .argsz = sizeof(group_status) };
	struct vfio_iommu_type1_info iommu_info = { .argsz = sizeof(iommu_info) };
	// VRP Firmware memory buffer
	struct vfio_iommu_type1_dma_map vrp_buffer_firmware = { .argsz = sizeof(vrp_buffer_firmware) };
	// VRP Data memory buffer
	struct vfio_iommu_type1_dma_map vrp_buffer_data = { .argsz = sizeof(vrp_buffer_data) };

	struct vfio_device_info device_info = { .argsz = sizeof(device_info) };

	/* Create a new container */
	*container = open(VFIO_CONTAINER, O_RDWR);

	if (ioctl(*container, VFIO_GET_API_VERSION) != VFIO_API_VERSION) {
		printf("Unknown API version\n");
		return -1;
	}

	if (!ioctl(*container, VFIO_CHECK_EXTENSION, VFIO_TYPE1_IOMMU)) {
		printf("Doesn't support the IOMMU driver we want\n");
		return -1;
	}

	/* Open the group */
	group = open(VRP_VFIO_GROUP, O_RDWR);

	/* Test the group is viable and available */
	ioctl(group, VFIO_GROUP_GET_STATUS, &group_status);

	if (!(group_status.flags & VFIO_GROUP_FLAGS_VIABLE)) {
		printf("Group is not viable (not all devices bound for vfio)\n");
		return -1;
	}

	/* Add the group to the container */
	ioctl(group, VFIO_GROUP_SET_CONTAINER, container);

	/* Enable the IOMMU model we want */
	ioctl(*container, VFIO_SET_IOMMU, VFIO_TYPE1_IOMMU);

	/* Get addition IOMMU info */
	ioctl(*container, VFIO_IOMMU_GET_INFO, &iommu_info);

	*size_to_map = getpagesize() * VRP_TRANSFER_SIZE_4KPAGES;

	// Get a buffer for the source of the VRP transfer and then map into the IOMMU

	vrp_buffer_firmware.vaddr = (u64)((uintptr_t)mmap(NULL, *size_to_map, PROT_READ | PROT_WRITE,
			     MAP_PRIVATE | MAP_ANONYMOUS, 0, 0));
    if ((void *)vrp_buffer_firmware.vaddr == MAP_FAILED) handle_error("mmap");

	vrp_buffer_firmware.size = *size_to_map;
	vrp_buffer_firmware.iova = vrp_buffer_firmware.vaddr;
	vrp_buffer_firmware.flags = VFIO_DMA_MAP_FLAG_READ | VFIO_DMA_MAP_FLAG_WRITE;

	l_rc = ioctl(*container, VFIO_IOMMU_MAP_DMA, &vrp_buffer_firmware);
	if(l_rc) {
		printf("Could not map VRP Firmware memory buffer\n");
		return -1;
	}
    *vrp_fw_mem_buff = vrp_buffer_firmware.vaddr;
    printf("[VFIO] VRP Firmware buffer: mapped %lld bytes at address 0x%lx\n", vrp_buffer_firmware.size, *vrp_fw_mem_buff);

	// Get a buffer for the destination of the VRP transfer and then map it into the IOMMU

	vrp_buffer_data.vaddr = (u64)((uintptr_t)mmap(NULL, *size_to_map, PROT_READ | PROT_WRITE,
			     MAP_PRIVATE | MAP_ANONYMOUS, 0, 0));
    if ((void *)vrp_buffer_data.vaddr == MAP_FAILED) handle_error("mmap");
	vrp_buffer_data.size = *size_to_map;
	vrp_buffer_data.iova = vrp_buffer_data.vaddr;
	vrp_buffer_data.flags = VFIO_DMA_MAP_FLAG_READ | VFIO_DMA_MAP_FLAG_WRITE;

	l_rc = ioctl(*container, VFIO_IOMMU_MAP_DMA, &vrp_buffer_data);
	if(l_rc) {
		printf("Could not map VRP Data memory buffer\n");
		return -1;
	}
    *vrp_data_mem_buff = vrp_buffer_data.vaddr;
    printf("[VFIO] VRP Data buffer: mapped %lld bytes at address 0x%lx\n", vrp_buffer_data.size, *vrp_data_mem_buff);

	/* Get a file descriptor for the device */
	device = ioctl(group, VFIO_GROUP_GET_DEVICE_FD, VRP_DEVICE);
	printf("=== VFIO device file descriptor %d ===\n", device);

	/* Test and setup the device */
	l_rc = ioctl(device, VFIO_DEVICE_GET_INFO, &device_info);

	if(l_rc) {
		printf("Could not get VFIO device\n");
		return -1;
	}

	printf("Device has %d region(s):\n", device_info.num_regions);

	struct vfio_region_info reg = { .argsz = sizeof(reg) };

	//uchar *base_vrp_reg;
	//uchar *base_rst_reg;
	//uchar *base_axi_reg;

    // VRP MAP
	reg.index = 0;
	l_rc = ioctl(device, VFIO_DEVICE_GET_REGION_INFO, &reg);

	if(l_rc) {
		printf("Couldn't get region %d info\n", reg.index);
		return -1;
	}

	printf("- Region %d: size=0x%llx offset=0x%llx flags=0x%x\n",
			reg.index,
			reg.size,
			reg.offset,
			reg.flags );

	*base_vrp_reg = (uint64_t)mmap(NULL, reg.size, PROT_READ | PROT_WRITE, MAP_SHARED,
			device, reg.offset);

	if ((void *)(*base_vrp_reg) != MAP_FAILED) {
		printf("Successful MMAP of VRP to address %p\n", (void *)*base_vrp_reg);
    } else {
		printf("Failed MMAP of VRP\n");
        return -1;
    }

    // VRP RESET GPIO
	reg.index = 1;
	l_rc = ioctl(device, VFIO_DEVICE_GET_REGION_INFO, &reg);

	if(l_rc) {
		printf("Couldn't get region %d info\n", reg.index);
		return -1;
	}

	printf("- Region %d: size=0x%llx offset=0x%llx flags=0x%x\n",
			reg.index,
			reg.size,
			reg.offset,
			reg.flags );

	*base_rst_reg = (uint64_t)mmap(NULL, reg.size, PROT_READ | PROT_WRITE, MAP_SHARED,
			device, reg.offset);

	if ((void *)*base_rst_reg != MAP_FAILED) {
		printf("Successful MMAP of RST_GPIO to address %p\n", (void *)*base_rst_reg);
    } else {
		printf("Failed MMAP of RST_GPIO\n");
        return -1;
    }

    // VRP AXI CACHE/PROT GPIO
	reg.index = 2;
	l_rc = ioctl(device, VFIO_DEVICE_GET_REGION_INFO, &reg);

	if(l_rc) {
		printf("Couldn't get region %d info\n", reg.index);
		return -1;
	}

	printf("- Region %d: size=0x%llx offset=0x%llx flags=0x%x\n",
			reg.index,
			reg.size,
			reg.offset,
			reg.flags );

	*base_axi_reg = (uint64_t)mmap(NULL, reg.size, PROT_READ | PROT_WRITE, MAP_SHARED,
			device, reg.offset);

	if ((void *)*base_axi_reg != MAP_FAILED) {
		printf("Successful MMAP of AXI_GPIO to address %p\n", (void *)*base_axi_reg);
    } else {
		printf("Failed MMAP of AXI_GPIO\n");
        return -1;
    }

    struct vfio_irq_set* vrp_irq_set = (struct vfio_irq_set*)malloc(sizeof(struct vfio_irq_set)+(sizeof(uint32_t)));
    struct vrp_irq_efd* vrp_efd_set = (struct vrp_irq_efd*)malloc(sizeof(struct vrp_irq_efd)+(sizeof(uint32_t)*device_info.num_irqs));
    vrp_efd_set->num_irqs = device_info.num_irqs;
    vrp_irq_set->argsz = sizeof(struct vfio_irq_set)+(sizeof(uint32_t));
    printf("Found %d VRP Interrups\n", device_info.num_irqs);
    // Map VRP IRQ to interrupt handler
    for (i = 0; i < device_info.num_irqs; i++) {
        struct vfio_irq_info irq = { .argsz = sizeof(irq) };
        irq.index = i;
        ioctl(device, VFIO_DEVICE_GET_IRQ_INFO, &irq);
        /* Setup IRQs... eventfds, VFIO_DEVICE_SET_IRQS */
        vrp_efd_set->irq_efd[i] = eventfd(0, EFD_NONBLOCK);
        if (vrp_efd_set->irq_efd[i] == (uint32_t)-1) l_rc = -1;
        *((uint32_t*)(vrp_irq_set->data)) = vrp_efd_set->irq_efd[i];
        vrp_irq_set->flags = VFIO_IRQ_SET_DATA_EVENTFD | VFIO_IRQ_SET_ACTION_TRIGGER;
        vrp_irq_set->index = i;
        vrp_irq_set->start = 0;
        vrp_irq_set->count = 1;

        ioctl(device, VFIO_DEVICE_SET_IRQS, vrp_irq_set);
    }

    *vrp_efd = vrp_efd_set;

    return l_rc;
}

int vrp_vfio_iommu_unmap(int* container, int* size_to_map, struct vrp_irq_efd** vrp_efd) {
    int l_rc = 0;
    struct vfio_iommu_type1_dma_unmap vrp_unmap_firmware = { .argsz = sizeof(vrp_unmap_firmware) };
	struct vfio_iommu_type1_dma_unmap vrp_unmap_data = { .argsz = sizeof(vrp_unmap_data) };

	vrp_unmap_firmware.size = *size_to_map;
	vrp_unmap_firmware.iova = 0;

	vrp_unmap_data.size = *size_to_map;
	vrp_unmap_data.iova = vrp_unmap_firmware.size;

    struct vfio_irq_set  vrp_irq_set;

    // UNMAP IRQ
    vrp_irq_set.start = 0;
    vrp_irq_set.count = 0;
    vrp_irq_set.flags = VFIO_IRQ_SET_DATA_NONE | VFIO_IRQ_SET_ACTION_TRIGGER;
    for (uint32_t i = 0; i < (*vrp_efd)->num_irqs; i++) {
        vrp_irq_set.index = i;
        ioctl(*container, VFIO_DEVICE_SET_IRQS, &vrp_irq_set);
        close((*vrp_efd)->irq_efd[i]);
    }
    free(*vrp_efd);
    *vrp_efd = NULL;

	l_rc = ioctl(*container, VFIO_IOMMU_UNMAP_DMA, &vrp_unmap_firmware);
	if(l_rc) {
		printf("Could not unmap VRP Firmware memory buffer\n");
		return -1;
	}

	l_rc = ioctl(*container, VFIO_IOMMU_UNMAP_DMA, &vrp_unmap_data);
	if(l_rc) {
		printf("Could not unmap VRP Data memory buffer\n");
		return -1;
	}

    return l_rc;
}

