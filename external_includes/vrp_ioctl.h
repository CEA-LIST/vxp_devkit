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
 * Creation Date : June, 2022
 * Description   : 
 **/
#ifndef __VRP_IOCTL_H__
#define __VRP_IOCTL_H__

#define VRP_IOCTL_TYPE              666

#define VRP_LOAD_FIRMWARE           _IOW(VRP_IOCTL_TYPE, 1, void *)
#define VRP_RUN_FIRMWARE            _IOW(VRP_IOCTL_TYPE, 2, void *)
#define VRP_STOP_FIRMWARE           _IOW(VRP_IOCTL_TYPE, 3, void *)
#define VRP_SET_SOLVER_ARGUMENTS    _IOW(VRP_IOCTL_TYPE, 4, void *)
#define VRP_GET_SOLVER_ARGUMENTS    _IOW(VRP_IOCTL_TYPE, 5, void *)
#define VRP_REGISTER_CALLER_APP     _IOW(VRP_IOCTL_TYPE, 6, void *)
#define VRP_UNREGISTER_CALLER_APP   _IOW(VRP_IOCTL_TYPE, 7, void *)
#define VRP_INIT_MEMORY             _IOW(VRP_IOCTL_TYPE, 8, void *)

typedef struct vrp_load_firmware_command { 
    char * firmware_file_path;          /**< Path to the firmware .bin file */
    size_t firmware_file_path_length;   /**< Length of the path string */
    size_t firmware_data_length;        /**< Size of the firmware so download */
} vrp_load_firmware_command_t;

typedef struct vrp_set_solver_data_cmd {
    uint64_t ndiag;
    double * matrix;
} vrp_set_solver_data_cmd_t;


typedef uint64_t vrp_count_t;
typedef uint64_t vrp_variable_size_t;
typedef uint64_t vrp_address_offset_t;
typedef enum vrp_argument_direction {
    VRP_SOLVER_ARGUMENT_IN,
    VRP_SOLVER_ARGUMENT_OUT,
    VRP_SOLVER_ARGUMENT_IN_OUT
} vrp_argument_direction_t;

/**
 * @brief 
 * 
 */
typedef struct vrp_solver_argument {
    vrp_argument_direction_t direction;
    vrp_variable_size_t size;
    vrp_variable_size_t alignment;
    void * address;
} vrp_solver_argument_t;

/**
 * @brief 
 * 
 */
typedef struct vrp_solver_argument_array{
    vrp_count_t nb_arguments;
    vrp_solver_argument_t arguments[];
} vrp_solver_argument_array_t;

#define sizeof_vrp_solver_arguments(x) (sizeof(vrp_count_t) + sizeof(vrp_solver_argument_t) * x)

typedef struct vrp_caller_application {
    int eventfd;
} vrp_caller_application_t;

/**
 * @brief 
 * 
 */
typedef struct vrp_init_memory_argument {
    vrp_variable_size_t user_buffer_size;
    void * user_buffer_address;
    vrp_count_t nb_write;
    uint64_t destination_start_address;
} vrp_init_memory_argument_t;


#endif /* __VRP_IOCTL_H__ */