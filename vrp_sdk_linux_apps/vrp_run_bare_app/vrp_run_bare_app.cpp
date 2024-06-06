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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "VRPOffload/vrp_driver_interface.hpp"

int main(int argc, char ** argv) {

    if ( argc != 3 ) {
        printf("Usage: %s <fw_directory> <fw_file_name>\n", argv[0]);
        exit(1);
    }

    char * l_firmware_dirname = argv[1];
    char * l_firmware_filename = argv[2];

    setenv("VRP_SOLVER_WRAPPERS_PATH", l_firmware_dirname, 1);
    setenv("VRP_OFFLOAD", "1", 1);

    printf("l_firmware_filename : %s\n", l_firmware_filename);
    printf("l_firmware_dirname : %s\n", l_firmware_dirname);

    if ( initialize_memory(1024*1024*1024, 0x800000000000) ) {
        printf("Fail initializing VRP memory.\n");
        exit(1);
    }

    if ( load_vrp_firmware(l_firmware_filename) != 0 ) {
        printf("Fail loading firmware in VRP memory.\n");
        exit(1);
    }

    start_vrp();

}