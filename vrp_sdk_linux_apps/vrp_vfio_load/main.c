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
 *  @file   main.c : Map the VRP in the SMMU and lunch VRP applications
 *  @author Riccardo ALIDORI
 *
******************************************************************************/
// NOTE:
// The vfio driver is used with the following commands prior to running this test application.
// Maybe it can be used from the device tree like UIO, but not sure yet.

// modprobe vfio_platform reset_required=0
// echo 20100000000.vrp > /sys/bus/platform/drivers/vrp_iommu/unbind
// echo vfio-platform > /sys/bus/platform/devices/20100000000.vrp/driver_override
// echo 20100000000.vrp > /sys/bus/platform/drivers_probe
// ******************************************************************************

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <libgen.h>
#include <dlfcn.h>
#include <signal.h>
#include <pthread.h>

#include "vrp_ioctl.h"
#include "vrp_vfio.h"
#include "vrp_init.h"
#include "vrp_irq.h"

#define MAX_LIMIT 50

static volatile int keepRunning = 1;

void intHandler(void) {
    printf("Catch CTRL-C, Free and Exit\n");
    keepRunning = 0;
}

void* vrp_irq_handler(void* p){
    struct vrp_irq_efd *vrp_efd = (struct vrp_irq_efd*)p;
    uint64_t u = 0;
    ssize_t s;
    int irq_id;
    uint64_t vrp_msi_msg;

    printf("Waiting for %d VRP IRQ\n", (*vrp_efd).num_irqs);
    while (keepRunning) {
        for (irq_id = 0; irq_id < (*vrp_efd).num_irqs; irq_id++) {
            s = read((*vrp_efd).irq_efd[irq_id], &u, sizeof(uint64_t));
            if (s == sizeof(uint64_t)) {

                printf("GIC : Got IRQ %d !\n",irq_id);
                /* Add your IRQ handler here */ //TODO pour le cher Jerome :)

                // Sending MSI to VRP with message 0xCEA + IRQ number
                vrp_msi_msg = 0xCEA + irq_id;
                printf("Sending MSI to VRP with message %d !\n",vrp_msi_msg);
                vrp_msi_it_send(vrp_msi_msg);

                u = 0;
                continue;
            }
        }
    }
    pthread_exit(EXIT_SUCCESS);
}

int vrp_submit_ioctl(unsigned long cmd, void * cmd_data) {
    int l_rc = 0;
    int l_vrp_fd = 0;

    l_vrp_fd = open("/dev/vrp0", O_RDWR);

    if ( l_vrp_fd == -1 ) {
        printf("Fail opening /dev/vrp0. %s\n", strerror(errno));
        return -ENODEV;
    }

    l_rc = ioctl(l_vrp_fd, cmd, cmd_data);
    if ( l_rc != 0 ) {
        printf("IOCTL %ld Failed!\n", cmd);
    }

    close(l_vrp_fd);
    return l_rc;
}

int burn_vrp(char * a_firmware_path, void * a_firmware_data, size_t a_firmware_size) {

    vrp_load_firmware_command_t l_load_firmware_cmd;
    l_load_firmware_cmd.firmware_path = a_firmware_path;
    l_load_firmware_cmd.firmware_data_address = a_firmware_data;
    l_load_firmware_cmd.firmware_data_length = a_firmware_size;

    return vrp_submit_ioctl(VRP_LOAD_FIRMWARE, &l_load_firmware_cmd);

}

int vrp_wait_for_solver_completion(uint64_t *vrp_fw_mem_buff) {
    uint64_t *shutdown_ptr;
	uint64_t l_solver_status;
    int timeout = 50;

    shutdown_ptr = vrp_fw_mem_buff;
	l_solver_status = *shutdown_ptr;

	while ( l_solver_status == VRP_SOLVER_STATUS_SUBMITED && timeout > 0 && keepRunning == 1) {
        usleep(500000); // 0.5sec
	    l_solver_status = *shutdown_ptr;
        //timeout--;
	}

    if (timeout == 0) {
        printf("[VRP IOMMU] Solver timeout!\n");
        return -1;
    }
	printf("[VRP IOMMU] l_solver_status %lld.\n", l_solver_status);

	return ( l_solver_status != VRP_SOLVER_STATUS_COMPLETED );
}

int vrp_dynamic_link(void* l_fw_data_copy, char * l_fw_path, int ram_size) {
    int l_rc = 0;
    int wdt = 0;
    char str_ram_uncached[50];
    char str_ram_cached[50];
    char str_ram_size[50];
    char str_path[100];
    char cmd[500];

    sprintf(str_ram_uncached,"%016llx", l_fw_data_copy);
    sprintf(str_ram_cached,  "%016llx", (((uint64_t)(l_fw_data_copy))+0x200000));
    sprintf(str_ram_size,  "%d", ram_size);
    sprintf(str_path,  "%s", l_fw_path);
    sprintf(cmd,  "%s %s %s %s %s",
            VRP_DYN_LINK_SCRIPT_PATH, // Script path
            str_ram_uncached,         // __ADDR_RAM_UNCACHED__
            str_ram_cached,           // __ADDR_RAM_CACHED__
            str_ram_size,             // __RAM_SIZE__ (MB)
            l_fw_path                 // Path of the file.x.a
            );

    printf("Call the dynamic linker script\n");
    printf("%s\n", cmd);

    //Call the dynamic linker script
    l_rc = system(cmd);

    if ( l_rc != 0 ) {
        printf("[VRP]Fail calling the dynamic linker script\n");
        return l_rc;
    }

    // Wait the _done file is created or the timeout (2mins) before continuing
    do {
        usleep(500000); // 0.5sec
        wdt++;
        if (wdt == 240) l_rc = 1;
    } while(access(VRP_DYN_LINK_DONE_PATH, F_OK) && (l_rc == 0));

    if ( l_rc != 0 ) {
        printf("[VRP]Fail Execution timeout dynamic Linking script\n");
    }

    return l_rc;
}

int load_vrp_firmware(void * vrp_fw_mem_buff, const char * a_firmware_path, int size) {
    int l_rc = 0;
    struct stat l_firmware_stat;
    int l_firmware_fd = -1;
    void * l_firmware_data = NULL;
    void * vrp_fw_handle;
    char * error;
    char * l_firmware_name;
    char buffer[512];
    int bytes_read;
    int off = 0;

    /* TODO restore this code
    if ( !getenv("VRP_SOLVER_WRAPPERS_PATH") ) {
        handle_error("VRP_SOLVER_WRAPPERS_PATH is not set.");
    }

    char * l_vrp_solver_wrappers_path = getenv("VRP_SOLVER_WRAPPERS_PATH");
    int l_firmware_path_length = ( strlen(a_firmware_path) + strlen(l_vrp_solver_wrappers_path) + 2 ); // +1 for additional '/' and +1 for leading \0
    char * l_firmware_path = (char *)malloc( sizeof(char)  *  l_firmware_path_length );

    if ( l_firmware_path == NULL ) {
        handle_error("Fail allocating memory to store firmware path");
    }

    snprintf(l_firmware_path, l_firmware_path_length, "%s/%s", l_vrp_solver_wrappers_path, a_firmware_path);
    */

    char * l_firmware_path = (char *)malloc( strlen(a_firmware_path) +1 );
    strcpy(l_firmware_path, a_firmware_path); // TODO Remove this afterwards

    if ( access(l_firmware_path, F_OK)) {
        handle_error("access");
    }

    // Dynamic Linking the VRP FW
    l_rc = vrp_dynamic_link(vrp_fw_mem_buff, a_firmware_path, ((size/1024)/1024));
    if (l_rc) {
        fprintf(stderr, "Error: Failed Dynamic linking the VRP FW\n");
        return -1;
    }

    l_firmware_fd  = open(VRP_FW_BINARY_PATH, O_RDONLY);

    if ( l_firmware_fd == -1 ) {
        handle_error("open");
    }

    // retrieve file size
    if ( fstat(l_firmware_fd, &l_firmware_stat) == -1 ) {
        handle_error("fstat");
    }

    // Copy the VRP binary in memory
    do {
        bytes_read = read(l_firmware_fd, buffer, 512);
        memcpy(vrp_fw_mem_buff + off, buffer, bytes_read);
        off += bytes_read;
    }
    while (bytes_read != 0);

    free(l_firmware_path);
    close(l_firmware_fd);

    return l_rc;
}

char * get_file_path() {
    char *str;
    str = malloc(sizeof(char)*MAX_LIMIT);

    printf("[VRP INIT] Charge the VRP FW automatically from the rootfs /usr/share/vrp/hello_world.x.a ? [y/n]:\n");
    fgets(str, MAX_LIMIT, stdin);
    // If yes
    if (!strcmp(str, "y\n")) {
        strcpy(str, "/usr/share/vrp/hello_world.x.a");
    // If no
    } else {
        printf("[VRP INIT] Insert VRP firmware absolute path:\n");

        if (fgets(str, MAX_LIMIT, stdin) == NULL) {
            printf("Error in getting string path\n");
            return "NULL";
        }
        // Remove trailing \n character
        str[strcspn(str, "\n")] = 0;
    }

    printf("return %s\n", str);
    return str;
}

int kick_vrp(uint64_t base_vrp_reg, uint64_t base_rst_reg, uint64_t base_axi_reg, uint64_t vrp_fw_mem_buff) {
    uint64_t *shutdown_ptr;
	printf("[VRP IOMMU] Setting shutdown pointer.\n");
    shutdown_ptr = vrp_fw_mem_buff;
    *shutdown_ptr = VRP_SOLVER_STATUS_SUBMITED;
    vrp_reset_platform(base_rst_reg);
    vrp_init_core(base_vrp_reg, base_rst_reg, base_axi_reg, vrp_fw_mem_buff);
    vrp_release_core(base_vrp_reg, base_rst_reg);
}

int main(int argc, char** argv) {
    char *path;
    const int ITER = 1;
    int container, size_to_map;
    uint64_t  vrp_fw_mem_buff;
    uint64_t  vrp_data_mem_buff;
    uint64_t  base_vrp_reg;
    uint64_t  base_rst_reg;
    uint64_t  base_axi_reg;
    //struct vfio_irq_set* vrp_irq_set = NULL;
    struct vrp_irq_efd *vrp_efd = NULL;
    fd_set rfds;
    struct timeval tv;
    /* Wait up to five seconds. */
    tv.tv_sec = 5;
    tv.tv_usec = 0; 

    // Catch Ctrl-C and stop VRP execution
    signal(SIGINT, intHandler);

    //
    // Map the VRP in the SMMU
    //
    printf("Mapping the VRP in the IOMMU...\n");
    if (vrp_vfio_iommu_map( &container,
                            &size_to_map,
                            &vrp_fw_mem_buff,
                            &vrp_data_mem_buff,
                            &base_vrp_reg,
                            &base_rst_reg,
                            &base_axi_reg,
                            &vrp_efd)) {
        perror("Failed to Map the VRP in the SMMU");
        return -1;
    }

    // Spawn a thread for polling on the VRP IRQs
    pthread_t tid;
    pthread_create(&tid, NULL, vrp_irq_handler, (void*)vrp_efd);

    //
    // Offload
    //

    //path = get_file_path();
    path = malloc(sizeof(char)*MAX_LIMIT);
    if (argc==2) {
        strcpy(path, argv[1]);
    } else {
        strcpy(path, "/usr/share/vrp/hello_world.x.a");
    }

    // Load VRP FW to Memory
    printf("[VRP INIT] Loading %s\n", path);
    if(load_vrp_firmware(vrp_fw_mem_buff, path, size_to_map)) {
        printf("Error in loading VRP firmware\n");
        goto exit;
    }

    // Kick VRP
    printf("[VRP INIT] Kick VRP\n");
    if (kick_vrp(base_vrp_reg, base_rst_reg, base_axi_reg, vrp_fw_mem_buff)) {
        printf("Error kicking the VRP\n");
        goto exit;
    }

    // Wait VRP completion
    printf("[VRP INIT] Wait for the VRP solver completion\n");
    if(vrp_wait_for_solver_completion(vrp_fw_mem_buff)) {
        printf("Error waiting for VRP completion\n");
        goto exit;
    }

exit:
    // Reset VRP Core before leaving
    vrp_reset_core(base_rst_reg);

    //
    // Unmap the VRP from the SMMU
    //
    free(path);
    printf("[VRP INIT] Unmapping the VRP from the IOMMU...\n");
    if (vrp_vfio_iommu_unmap(&container, &size_to_map, &vrp_efd)) {
        perror("Failed to Unmap the VRP in the SMMU");
        return -1;
    }

    // We signal the IRQ thread to stop working
    keepRunning = 0;

    printf("[VRP INIT] done.\n");
    return 0;
}
