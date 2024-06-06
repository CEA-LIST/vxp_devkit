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
#ifndef __VRP_IOCTL_H__
#define __VRP_IOCTL_H__

#define VRP_IOCTL_TYPE              666

#define VRP_LOAD_FIRMWARE           _IOW(VRP_IOCTL_TYPE, 1, void *)
#define VRP_SET_SOLVER_ARGUMENTS    _IOW(VRP_IOCTL_TYPE, 2, void *)
#define VRP_GET_SOLVER_ARGUMENTS    _IOW(VRP_IOCTL_TYPE, 3, void *)

typedef struct vrp_load_firmware_command {
    char * firmware_path;
    void * firmware_data_address;
    size_t firmware_data_length;
} vrp_load_firmware_command_t;

typedef struct vrp_set_solver_data_cmd {
    uint64_t ndiag;
    double * matrix;
} vrp_set_solver_data_cmd_t;

typedef uint64_t vrp_count_t;
typedef uint64_t vrp_variable_size_t;
typedef uint64_t vrp_address_offset_t;

/**
 * @brief 
 * 
 */
typedef struct vrp_solver_argument {
    vrp_variable_size_t size;
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

#define handle_error(msg) \
    do { perror(msg); exit(EXIT_FAILURE); } while (0)

#endif /* __VRP_IOCTL_H__ */
