/**
* Copyright 2022 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 *  @file        vrp_solver_qmr_firmware.c
 *  @author      Jerome Fereyre
 */

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <stdint.h>

#include "VRPOffload/vrp_Matrix_serializer.hpp"
#include "VPSDK/VBLASConfig.hpp"
#include "VPSolvers.hpp"

#include "alignment.h"
#include "VRPSDK/perfcounters/cpu.h"

void * __dso_handle = NULL;

using namespace VPFloatPackage;

int qmr_wrapper(){
    std::ostringstream strCout;

    int l_matrix_format_invalid = 0;

    uint64_t * l_vrp_solver_status_ptr = (uint64_t *)VRP_DATA_ADDRESS; 
    uint64_t l_param_address = VRP_DATA_ADDRESS + sizeof(uint64_t);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    int precision = *(int *)(l_param_address);
    l_param_address += sizeof(int);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    int transpose = *(int *)(l_param_address);
    l_param_address += sizeof(int);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    int n = *(int *)(l_param_address);
    l_param_address += sizeof(int);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    double tolerance = *(double *)(l_param_address);
    l_param_address += sizeof(double);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    double * X = (double *)(l_param_address);
    l_param_address += sizeof(double) * n;
    ALIGN_ADDRESS_64Bytes(l_param_address);

    matrix_t A = Offloading::VRP_Matrix_serializer::unserialize(&l_param_address);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    matrix_t At = Offloading::VRP_Matrix_serializer::unserialize(&l_param_address);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    double * B = (double *)(l_param_address);
    l_param_address += sizeof(double) * n;
    ALIGN_ADDRESS_64Bytes(l_param_address);
 
    int16_t exponent_size = *(int16_t *)(l_param_address);
    l_param_address += sizeof(int16_t);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    int32_t stride_size = *(int32_t *)(l_param_address);
    l_param_address += sizeof(int32_t);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    int32_t * iteration_count = (int32_t *)(l_param_address);
    l_param_address += sizeof(int32_t);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    uint64_t log_buffer_size = *(uint64_t *)(l_param_address);
    l_param_address += sizeof(uint64_t);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    char * log_buffer = (char *)(l_param_address);
    l_param_address += sizeof(char) * (log_buffer_size);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    if ( log_buffer_size >= 0 ) {
      printf("log_buf_size is %ld. Redirect cout to oStringStream.\n", log_buffer_size);
      std::cout.rdbuf( strCout.rdbuf() );
    }

    VBLAS::VBLAS_setConfig((VBLAS::VBLASConfig *)l_param_address);
    l_param_address += sizeof(VBLAS::VBLASConfig);
    ALIGN_ADDRESS_64Bytes(l_param_address);

    std::cout<<"1.QMR vanille , ";

    l_matrix_format_invalid = displayMatrixCharacteristics(A);

    std::cout << std::endl;

    if ( l_matrix_format_invalid ) {
      *l_vrp_solver_status_ptr = (uint64_t)0x0000000000000002;
      return 1;
    }

    uint64_t cy=get_cycles();
    uint64_t l_start_nb_instructions = cpu_instructions();

    int32_t l_nb_iteration = VPFloatPackage::Solver::qmr(precision, transpose, n, X, A, At, B, tolerance, exponent_size, stride_size, log_buffer, log_buffer_size);

    uint64_t nbcycles=get_cycles()-cy;
    uint64_t l_nb_instructions = cpu_instructions() - l_start_nb_instructions;

    std::cout << precision << " "<< l_nb_iteration << " " << nbcycles<<" " <<  double(l_nb_instructions) / nbcycles << "\n";

    *iteration_count=l_nb_iteration;

    if ( log_buffer_size >= 0 ) {
      strncpy(log_buffer, strCout.str().c_str(), std::min(strCout.str().length(), log_buffer_size));
    }

    *l_vrp_solver_status_ptr = (uint64_t)0x0000000000000003;

    return 0;

}

int main(int argc, char *argv[])
{
  printf("Calling qmr_wrapper.\n");
  
  int l_rc = qmr_wrapper();

  exit(l_rc);
}
