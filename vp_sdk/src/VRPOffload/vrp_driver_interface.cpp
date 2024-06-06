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

#include <string.h>
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
#include <iostream>
#include <time.h>
#include <signal.h>
#include <sys/eventfd.h>


#include "VRPOffload/vrp_driver_interface.hpp"
#include "VRPOffload/vrp_offloading.hpp"
#include "vrp_ioctl.h"

using namespace VPFloatPackage::Offloading;

struct sigaction g_previous_sigint_action;


int vrp_submit_ioctl(unsigned long cmd, void * cmd_data) {
    int l_rc = 0;
    int l_vrp_fd = 0;

    l_vrp_fd = open("/dev/vrp0", O_RDWR);

    if ( l_vrp_fd == -1 ) {
        printf("Fail opening /dev/vrp0. %s\n", strerror(errno));
        return -ENODEV;
    }
    if ( cmd_data == NULL ) {
        l_rc = ioctl(l_vrp_fd, cmd);
    } else {
        l_rc = ioctl(l_vrp_fd, cmd, cmd_data);
    }

    if ( l_rc != 0 ) {
        printf("IOCTL %ld Failed!\n", cmd);
    }

    close(l_vrp_fd);
    return l_rc;
}

void vrp_sigint_handler(int a_signal) {
    printf("Calling %s for signal %s.\n", __func__, strsignal(a_signal));
    vrp_submit_ioctl(VRP_STOP_FIRMWARE, NULL);

    if (    g_previous_sigint_action.sa_handler != SIG_IGN && 
            g_previous_sigint_action.sa_handler != SIG_DFL && 
            g_previous_sigint_action.sa_handler != NULL ) {
        printf("Call previously registered SIGTERM handler.\n");
        g_previous_sigint_action.sa_handler(a_signal);
    }
}

int burn_vrp(char * a_firmware_file_path, size_t a_firmware_size) {
    
    vrp_load_firmware_command_t l_load_firmware_cmd;
    l_load_firmware_cmd.firmware_file_path_length = strlen(a_firmware_file_path); 
    l_load_firmware_cmd.firmware_file_path = a_firmware_file_path; 
    l_load_firmware_cmd.firmware_data_length = a_firmware_size;

    return vrp_submit_ioctl(VRP_LOAD_FIRMWARE, &l_load_firmware_cmd);
    
}

int start_vrp() {
    
    int l_vrp_eventfd = eventfd(0, 0);
    uint64_t l_vrp_event_data;
    ssize_t l_vrp_event_data_size;
    vrp_caller_application_t * l_caller_application = (vrp_caller_application_t *) malloc(sizeof(vrp_caller_application_t));

    if( l_caller_application == NULL ) {
        handle_error("Fail allocating memory for caller application structure.\n");
    }

    if ( l_vrp_eventfd == -1 ) {
        free(l_caller_application);
        handle_error("Fail creating new eventfd for VRP driver notification.\n");
    }

    l_caller_application->eventfd = l_vrp_eventfd;

    int l_rc = vrp_submit_ioctl(VRP_REGISTER_CALLER_APP, (void *)l_caller_application);

    if ( l_rc != 0 ) {
        free(l_caller_application);
        handle_error("Fail registering application into VRP driver.\n");
    }

    l_rc = vrp_submit_ioctl(VRP_RUN_FIRMWARE, NULL);

    l_vrp_event_data_size = read(l_vrp_eventfd, &l_vrp_event_data, sizeof(uint64_t));

    if (l_vrp_event_data_size != sizeof(uint64_t)) {
        free(l_caller_application);
        handle_error("Fail reading data from eventfd for VRP driver notification.\n");
    }

    l_rc = vrp_submit_ioctl(VRP_UNREGISTER_CALLER_APP, l_caller_application);

    if ( l_rc != 0 ) {
        free(l_caller_application);
        handle_error("Fail unregistering application from VRP driver.\n");
    }

    free(l_caller_application);

    return l_rc;
}

int initialize_memory(size_t a_size, uint64_t a_destination_address) {
    size_t l_init_mem_size = ( a_size / 1024 ) * sizeof(char);
    vrp_init_memory_argument_t l_init_mem_cmd;

    // Allocate 100Mo buffer
    void * l_init_mem = malloc(l_init_mem_size);

    // Check that memory allocation succeed
    if ( l_init_mem == NULL ) {
        handle_error("Fail allocating init mem buffer.\n");
    }

    // Set memory to 0
    memset(l_init_mem, 0x0, l_init_mem_size);

    // Build the initialize memory command
    l_init_mem_cmd.destination_start_address = a_destination_address;
    l_init_mem_cmd.nb_write = 1024;
    l_init_mem_cmd.user_buffer_address = l_init_mem;
    l_init_mem_cmd.user_buffer_size = l_init_mem_size;

    int l_rc = vrp_submit_ioctl(VRP_INIT_MEMORY, &l_init_mem_cmd);

    // Release buffer used by memory initialization procedure
    free(l_init_mem);

    return l_rc;
}

int load_vrp_firmware(const char * a_firmware_path) {
    int l_rc = 0;
    struct stat l_firmware_stat;
    int l_firmware_fd = -1;
    void * l_firmware_data = NULL;

    if ( !getenv("VRP_SOLVER_WRAPPERS_PATH") ) {
        handle_error("VRP_SOLVER_WRAPPERS_PATH is not set.");
    }

    char * l_vrp_solver_wrappers_path = getenv("VRP_SOLVER_WRAPPERS_PATH");
    int l_firmware_path_length = ( strlen(a_firmware_path) + strlen(l_vrp_solver_wrappers_path) + 2 ); // +1 for additional '/' and +1 for leading \0
    char * l_firmware_path = (char *)malloc( sizeof(char)  *  l_firmware_path_length);
    memset(l_firmware_path, '\0', l_firmware_path_length);

    if ( l_firmware_path == NULL ) {
        handle_error("Fail allocating memory to store firmware path");
    }

    snprintf(l_firmware_path, l_firmware_path_length, "%s/%s", l_vrp_solver_wrappers_path, a_firmware_path);

    if ( access(l_firmware_path, F_OK)) {
        handle_error("access");
    }

    l_firmware_fd  = open(l_firmware_path, O_RDONLY);

    if ( l_firmware_fd == -1 ) {
        handle_error("open");
    }

    // retrieve file size
    if ( fstat(l_firmware_fd, &l_firmware_stat) == -1 ) {
        handle_error("fstat");
    }

    printf("loading firmware %s in VRP\n", l_firmware_path);

    l_rc = burn_vrp(l_firmware_path, l_firmware_stat.st_size);

    free(l_firmware_path);

    return l_rc;
}

int write_solver_arguments(vrp_solver_argument_array_t * a_solver_arguments) {
    return vrp_submit_ioctl(VRP_SET_SOLVER_ARGUMENTS, (void *)a_solver_arguments);
}

int read_solver_arguments(vrp_solver_argument_array_t * a_solver_arguments) {
    return vrp_submit_ioctl(VRP_GET_SOLVER_ARGUMENTS, (void *)a_solver_arguments);
}


namespace VPFloatPackage::Offloading {

    int call_solver(const VRPArgumentArray & a_argument_array, const char * a_solver_bin_path) {
        int l_rc = 0;
        struct timespec l_timespec_start, l_timespec_stop;
        uint64_t l_write_duration, l_read_duration, l_solver_duration, l_firmware_transfert_duration;
        struct sigaction l_vrp_sigterm_action;
        vrp_solver_argument_array_t * l_marshalled_argument_array = NULL;

        clock_gettime(CLOCK_MONOTONIC, &l_timespec_start);
        l_marshalled_argument_array = marshall(a_argument_array);
        l_rc = write_solver_arguments(l_marshalled_argument_array);

        if ( l_rc != 0 ) {
            printf("Fail writing solver argument into VRP.\n");
            return l_rc;
        }
        clock_gettime(CLOCK_MONOTONIC, &l_timespec_stop);

        l_write_duration = ( ( ( l_timespec_stop.tv_sec - l_timespec_start.tv_sec ) * 1e9 ) + ( l_timespec_stop.tv_nsec - l_timespec_start.tv_nsec ) );

        // This is a blocking call. Waiting for solver completion.. or error.
        clock_gettime(CLOCK_MONOTONIC, &l_timespec_start);
        l_rc = load_vrp_firmware(a_solver_bin_path);

        if ( l_rc != 0 ) {
            printf("Fail loading VRP firmware.\n");
            return l_rc;
        } 
        clock_gettime(CLOCK_MONOTONIC, &l_timespec_stop);

        l_firmware_transfert_duration =  ( ( ( l_timespec_stop.tv_sec - l_timespec_start.tv_sec ) * 1e9 ) + ( l_timespec_stop.tv_nsec - l_timespec_start.tv_nsec ) );

        // Save previous SIGTERM handler
        sigaction (SIGINT, NULL, &g_previous_sigint_action);
        
        // Register new SIGTERM handler
        l_vrp_sigterm_action.sa_handler = vrp_sigint_handler;
        l_vrp_sigterm_action.sa_flags = 0;
        sigemptyset(&l_vrp_sigterm_action.sa_mask);
        sigaction (SIGINT, &l_vrp_sigterm_action, NULL);

        clock_gettime(CLOCK_MONOTONIC, &l_timespec_start);
        l_rc = start_vrp();

        if ( l_rc != 0 ) {
            printf("Fail starting VRP firmware.\n");
            return l_rc;
        }
        clock_gettime(CLOCK_MONOTONIC, &l_timespec_stop);

        l_solver_duration = ( ( ( l_timespec_stop.tv_sec - l_timespec_start.tv_sec ) * 1e9 ) + ( l_timespec_stop.tv_nsec - l_timespec_start.tv_nsec ) );

        // Restore previsou SIGTERM handler
        sigaction (SIGINT, &g_previous_sigint_action, NULL);

        clock_gettime(CLOCK_MONOTONIC, &l_timespec_start);
        l_rc = read_solver_arguments(l_marshalled_argument_array);

        if ( l_rc != 0 ) {
            printf("Fail reading solver argument into VRP.\n");
            return l_rc;
        } 
        clock_gettime(CLOCK_MONOTONIC, &l_timespec_stop);

        free(l_marshalled_argument_array);

        l_read_duration = ( ( ( l_timespec_stop.tv_sec - l_timespec_start.tv_sec ) * 1e9 ) + ( l_timespec_stop.tv_nsec - l_timespec_start.tv_nsec ) );

        printf("firmware transfert duration : %ldns\n", l_firmware_transfert_duration);
        printf("solver duration             : %ldns\n", l_solver_duration);
        printf("write solver arg duration   : %ldns\n", l_write_duration);
        printf("read  solver arg duration   : %ldns\n", l_read_duration);

        return l_rc;
    }

}