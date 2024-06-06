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

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <time.h>
#include "bicg_kernel.hpp"
#include "VRPOffload/vrp_offloading.hpp"
#include "VPSDK/VBLASConfig.hpp"

using namespace VPFloatPackage::Offloading;

int bicg_with_vrp_offload(int precision, int transpose, int n, double * X, matrix_t A, matrix_t At, double * B, double tolerance, uint16_t exponent_size, int32_t stride_size, char * log_buffer, uint64_t log_buffer_size){
    int l_rc = 0;
    int l_iteration_count = 0;

    // Build solver argument array
    VRPArgumentArray l_argument_array;

    l_argument_array.addArgument(&precision, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&transpose, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&n, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&tolerance, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(X, n * sizeof(double), VRP_SOLVER_ARGUMENT_IN_OUT);
    l_argument_array.addArgument(A, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(At, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(B, n * sizeof(double), VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&exponent_size, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&stride_size, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&l_iteration_count, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(&log_buffer_size, VRP_SOLVER_ARGUMENT_IN);
    l_argument_array.addArgument(log_buffer, log_buffer_size * sizeof(char), VRP_SOLVER_ARGUMENT_IN_OUT);
    l_argument_array.addArgument(VBLAS::VBLAS_getConfig(), VRP_SOLVER_ARGUMENT_IN);

    l_rc = call_solver(l_argument_array, "vrp_solver_bicg.x.bin");

    if ( l_rc == 0 ) {
        return l_iteration_count;
    }

    return l_rc;
}

namespace VPFloatPackage::Solver {

    int bicg(int precision, int transpose, int n, double * x, matrix_t A, matrix_t At, double * b, double tolerance, uint16_t exponent_size, int32_t stride_size, char * log_buffer, uint64_t log_buffer_size) {
        char * l_vrp_offload = getenv(VRP_OFFLAD_ENVIRONMENT_VAR_NAME);

        if (l_vrp_offload != NULL && atoi(l_vrp_offload) != 0 ) {
            return  bicg_with_vrp_offload(precision, transpose, n, x, A, At, b, tolerance, exponent_size, stride_size, log_buffer, log_buffer_size);
        } else {
            std::streambuf *cout_backup_buf;
            std::ostringstream strCout;
            if ( log_buffer_size > 0 ) {
                printf("log_buf_size is %ld. Redirect cout to oStringStream.\n", log_buffer_size);
                cout_backup_buf = std::cout.rdbuf();
                std::cout.rdbuf( strCout.rdbuf() );
            }

            VPFloatArray Xv(x, n);
            VPFloatArray Bv(b, n);

            struct timespec l_timespec_start, l_timespec_stop;
            uint64_t l_solver_duration;

            clock_gettime(CLOCK_MONOTONIC, &l_timespec_start);
            int l_iteration_count =  bicg_vp(precision, transpose, n, Xv, A, At, Bv, tolerance, exponent_size, stride_size);
            clock_gettime(CLOCK_MONOTONIC, &l_timespec_stop);

            l_solver_duration = ( ( ( l_timespec_stop.tv_sec - l_timespec_start.tv_sec ) * 1e9 ) + ( l_timespec_stop.tv_nsec - l_timespec_start.tv_nsec ) );

            std::cout << "solver duration             : " << l_solver_duration << "ns" << std::endl;
            std::cout << precision << " "<< l_iteration_count << std::endl;

            if ( log_buffer_size >= 0 ) {
                strncpy(log_buffer, strCout.str().c_str(), std::min(strCout.str().length(), log_buffer_size));
                std::cout.rdbuf(cout_backup_buf);
            }

            return l_iteration_count;
        }
    }

}
