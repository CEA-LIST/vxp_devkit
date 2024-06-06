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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "VPSolvers.hpp"
#include "OSKIHelper.hpp"
#include "Preconditionners.hpp"
#include "MTXUtil/crs.h"
#include "VPSDK/VBLASConfig.hpp"
#include "VRPOffload/vrp_offloading.hpp"
#include "Matrix/DENSE.h"

using namespace VPFloatPackage::Solver;

void toUpper(char * a_string) {
    int i;
    for(i=0; a_string[i]!='\0'; i++) {
        if(a_string[i]>='a' && a_string[i]<='z') {
            a_string[i] = a_string[i] - 32;
        }
    }
}

void usage(int a_rc) {
    printf("-h : print help message.\n");
    printf("-a <lda_value>                          : padded size of matrice lines. Use for cache prefetching (default:0 => automatic LDA tunning)\n");
    printf("-b <block_size>                         : size for block in BCSR format\n");
    printf("-c                                      : enable hardware prefetching\n");
    printf("-e <exponent_size>                      : size of exponent for VPfloat number used during solver computation.(default: 10)\n");
    printf("-k <kernel_name>                        : name of the kernel to call (BICG, CG, PRECOND_CG, QMR).\n");
    printf("-m <matrix_path>                        : path to the matrix to which the selected solver will be applied.\n");
    printf("-o                                      : request solver offloading on VRP accelerator\n");
    printf("-p <precision>                          : precision used during solver computation. (default: 512)\n");
    printf("-s                                      : flag used to specify kernel work wirth sparse or dense data structure.\n");
    printf("-t <tolerance in scientific notation>   : tolerance used by solver to determine end of iteration.(default: 1e-8)\n");
    printf("-l <log buffer size in byte>            : size of the buffer given to VRP to store solver output traces.\n");
    printf("-y                                      : Request transposed version of algorithm to run.\n");
    exit(a_rc);
}

void initB(double * B, int n) {
    for ( int index = 0 ; index < n ; index++ ) {
        B[index] = 1.0;
    }
}

int main(int argc, char ** argv) {
    
    VPFloatPackage::VBLAS::VBLASConfig * l_vblas_config = VPFloatPackage::VBLAS::VBLAS_getConfig();

    int l_rc = 0;
    int l_opt;
    int l_precision = 512;
    int l_transpose = 0;
    char * l_matrix_file_path = NULL;
    char * l_B_matrix_file_path = NULL;
    char * l_solver_name = NULL;
    uint64_t l_log_buffer_size = 0;
    char * l_log_buffer = NULL;
    bool l_sparse_flag = false;
    double l_tolerance = 1e-8;
    uint16_t l_exponent_size = 10;
    uint16_t l_stride_size = 1;
    int l_lda = 0;
    bool l_bcsr_sparse_block_size = 0;
    matrix_t l_B_matrix_loaded_from_file = NULL;
    oski_matrix_wrapper_t l_oski_B_input_matrix;
    double l_jacobi_shifter = 0.0;

    // By default deactivate prefetcher
    l_vblas_config->enable_prefetcher = 0;

    while((l_opt = getopt(argc, argv, "a:B:b:ce:hj:k:l:m:on:p:r:st:y")) != -1 ) {
        switch(l_opt) {
            case 'a':
                l_lda = atoi(optarg);
                break;
            case 'B':
                l_B_matrix_file_path = (char *)malloc( ( strlen(optarg) + 2 ) * sizeof(char) );
                snprintf(l_B_matrix_file_path, strlen(optarg) + 1 , "%s", optarg);                
                break;                
            case 'b':
                l_bcsr_sparse_block_size = atoi(optarg);
                break;
            case 'c':
                l_vblas_config->enable_prefetcher = 1;
                break;                
            case 'e':
                sscanf(optarg, "%hd", &l_exponent_size);
                break;
            case 'h':
                usage(0);
                break;
            case 'j':
                sscanf(optarg, "%le", &l_jacobi_shifter);
                break;                 
            case 'o':
                setenv(VRP_OFFLAD_ENVIRONMENT_VAR_NAME, "1", 1);
                break;
            case 'p':
                l_precision = atoi(optarg);
                break;
            case 'm':
                l_matrix_file_path = (char *)malloc( ( strlen(optarg) + 2 ) * sizeof(char) );
                snprintf(l_matrix_file_path, strlen(optarg) + 1 , "%s", optarg);
                break;
            case 'n':
                l_vblas_config->nb_threads = atoi(optarg);
                break;                 
            case 'k':
                l_solver_name = (char *)malloc( ( strlen(optarg) + 2 ) * sizeof(char) );
                snprintf(l_solver_name, strlen(optarg) + 1 , "%s", optarg);
                toUpper(l_solver_name);
                break;
            case 'l':
                sscanf(optarg, "%ld", &l_log_buffer_size);
                break; 
            case 'r':
                l_vblas_config->nb_rows_per_thread = atoi(optarg);
                break;                
            case 's':
                l_sparse_flag = true;
                break;                
            case 't':
                sscanf(optarg, "%le", &l_tolerance);
                break;  
            case 'y':
                l_transpose = 1;
                break;                              
            default:
                usage(-1);
                break;                
        }
    }

    if ( l_matrix_file_path == NULL ) {
        printf("-m option is mandatory.\n Please specify a matrix to load.\n");
        exit(1);
    }

    printf("===== Tolerance set to %le.\n", l_tolerance);
    printf("===== Exponent size set to %hd.\n", l_exponent_size);

    oski_matrix_wrapper_t l_oski_sparse_input_matrix = VPFloatPackage::OSKIHelper::loadFromFile(l_matrix_file_path);

    if ( l_oski_sparse_input_matrix.oski_matrix.real_matrix == NULL && l_oski_sparse_input_matrix.oski_matrix.complex_matrix == NULL) {
        printf("Fail parsing matrix.\n");
        exit(1);
    }
    matrix_t l_sparse_input_matrix = VPFloatPackage::OSKIHelper::toMatrix(l_oski_sparse_input_matrix);

    oski_matrix_wrapper_t l_oski_sparse_input_matrix_transposed = VPFloatPackage::OSKIHelper::transpose(l_oski_sparse_input_matrix);
    matrix_t l_sparse_input_matrix_transposed = VPFloatPackage::OSKIHelper::toMatrix(l_oski_sparse_input_matrix_transposed);

    if ( l_log_buffer_size > 0 ) {
        l_log_buffer = (char *)malloc(sizeof(char) * l_log_buffer_size);
        memset(l_log_buffer, 0, sizeof(char) *l_log_buffer_size);
    }

    double * X = (double *)malloc(sizeof(double) * l_sparse_input_matrix->n);
    memset(X, 0, sizeof(double) * l_sparse_input_matrix->n);

    /* Using the right B vector if specified */
    double * B = NULL;
    if ( l_B_matrix_file_path == NULL ) {
        std::cout << "Use automatically generated B vector." << std::endl;
        B = (double *)malloc(sizeof(double) * l_sparse_input_matrix->n);
        memset(B, 0, sizeof(double) * l_sparse_input_matrix->n);

        initB(B, l_sparse_input_matrix->n);
    } else {
        std::cout << "Use B vector loaded from file." << std::endl;
        l_oski_B_input_matrix = VPFloatPackage::OSKIHelper::loadFromFile(l_B_matrix_file_path);
        l_B_matrix_loaded_from_file = VPFloatPackage::OSKIHelper::toDense(l_oski_B_input_matrix, false, 0);
        B = ((dmatDENSE_t)l_B_matrix_loaded_from_file->matrix->repr)->val;
    }
    
    if  ( l_sparse_flag ) {
        /*
         * SPARSE version of solvers
         */
        if (strcmp(l_solver_name, "BICG") == 0 || strcmp(l_solver_name, "BICGSTAB") == 0 || strcmp(l_solver_name, "PRECOND_BICG") == 0 ) {          
            if ( l_bcsr_sparse_block_size != 0) {
                matrix_t l_bcsr_input_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                matrix_t l_bcsr_input_matrix_transposed = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix_transposed, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);

                if (strcmp(l_solver_name, "BICG") == 0) {
                    l_rc = bicg(l_precision,
                                l_transpose,
                                l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n,  
                                X, 
                                l_transpose == 1 ? l_bcsr_input_matrix_transposed : l_bcsr_input_matrix, 
                                l_transpose == 1 ? l_bcsr_input_matrix : l_bcsr_input_matrix_transposed, 
                                B, 
                                l_tolerance, 
                                l_exponent_size, 
                                l_stride_size, 
                                l_log_buffer, 
                                l_log_buffer_size);
                } else if (strcmp(l_solver_name, "BICGSTAB") == 0) {
                    l_rc = bicgstab(l_precision, 
                                    l_transpose,
                                    l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                                    X, 
                                    l_transpose == 1 ? l_bcsr_input_matrix_transposed : l_bcsr_input_matrix, 
                                    l_transpose == 1 ? l_bcsr_input_matrix : l_bcsr_input_matrix_transposed, 
                                    B, 
                                    l_tolerance, 
                                    l_exponent_size, 
                                    l_stride_size, 
                                    l_log_buffer, 
                                    l_log_buffer_size);
                } else {
                    if ( l_transpose == 1 ) {
                        matrix_t l_bcsr_input_matrix_transposed = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix_transposed, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                        matrix_t l_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);

                        l_rc = precond_bicg(l_precision,
                                            l_transpose,
                                            l_sparse_input_matrix_transposed->n,  
                                            X, 
                                            l_bcsr_input_matrix_transposed, 
                                            l_bcsr_input_matrix, 
                                            l_iM_transposed,
                                            B, 
                                            l_tolerance, 
                                            l_exponent_size, 
                                            l_stride_size, 
                                            l_log_buffer, 
                                            l_log_buffer_size);
                    } else {
                        matrix_t l_bcsr_input_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                        matrix_t l_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);
                        
                        l_rc = precond_bicg(l_precision,
                                            l_transpose,
                                            l_sparse_input_matrix->n,  
                                            X, 
                                            l_bcsr_input_matrix, 
                                            l_bcsr_input_matrix_transposed, 
                                            l_iM,
                                            B, 
                                            l_tolerance, 
                                            l_exponent_size, 
                                            l_stride_size, 
                                            l_log_buffer, 
                                            l_log_buffer_size);
                    }
                }

            } else {                
                if (strcmp(l_solver_name, "BICG") == 0) {
                    l_rc = bicg(l_precision, 
                                l_transpose,
                                l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                                X, 
                                l_transpose == 1 ? l_sparse_input_matrix_transposed : l_sparse_input_matrix, 
                                l_transpose == 1 ? l_sparse_input_matrix : l_sparse_input_matrix_transposed, 
                                B, 
                                l_tolerance, 
                                l_exponent_size, 
                                l_stride_size, 
                                l_log_buffer, 
                                l_log_buffer_size);
                } else if (strcmp(l_solver_name, "BICGSTAB") == 0) {
                    l_rc = bicgstab(l_precision,
                                l_transpose,
                                l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                                X, 
                                l_transpose == 1 ? l_sparse_input_matrix_transposed : l_sparse_input_matrix, 
                                l_transpose == 1 ? l_sparse_input_matrix : l_sparse_input_matrix_transposed,
                                B, 
                                l_tolerance, 
                                l_exponent_size, 
                                l_stride_size, 
                                l_log_buffer, 
                                l_log_buffer_size);
                } else {
                    if ( l_transpose == 1 ) {
                        matrix_t l_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);

                        l_rc = precond_bicg(l_precision, 
                                            l_transpose,
                                            l_sparse_input_matrix_transposed->n, 
                                            X, 
                                            l_sparse_input_matrix_transposed, 
                                            l_sparse_input_matrix, 
                                            l_iM_transposed,
                                            B, 
                                            l_tolerance, 
                                            l_exponent_size, 
                                            l_stride_size, 
                                            l_log_buffer, 
                                            l_log_buffer_size);
                    } else {
                        matrix_t l_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);

                        l_rc = precond_bicg(l_precision, 
                                            l_transpose,
                                            l_sparse_input_matrix->n, 
                                            X, 
                                            l_sparse_input_matrix, 
                                            l_sparse_input_matrix_transposed, 
                                            l_iM,
                                            B, 
                                            l_tolerance, 
                                            l_exponent_size, 
                                            l_stride_size, 
                                            l_log_buffer, 
                                            l_log_buffer_size);
                    }
                }
            }
        } else if (strcmp(l_solver_name, "CG") == 0) {
            
            if ( l_bcsr_sparse_block_size != 0) {
                if ( l_transpose == 1 ) {
                    matrix_t l_bcsr_input_matrix_transposed = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix_transposed, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                    l_rc = cg(  l_precision,
                            l_transpose, 
                            l_sparse_input_matrix_transposed->n, 
                            X, 
                            l_bcsr_input_matrix_transposed,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size,
                            l_log_buffer, 
                            l_log_buffer_size);
                } else {
                    matrix_t l_bcsr_input_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                    l_rc = cg(  l_precision,
                            l_transpose, 
                            l_sparse_input_matrix->n, 
                            X, 
                            l_bcsr_input_matrix,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size,
                            l_log_buffer, 
                            l_log_buffer_size);
                }
            } else {
                l_rc = cg(  l_precision, 
                            l_transpose,
                            l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                            X, 
                            l_transpose == 1 ? l_sparse_input_matrix_transposed : l_sparse_input_matrix,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size,
                            l_log_buffer, 
                            l_log_buffer_size);
            }
        } else if(strcmp(l_solver_name, "PRECOND_CG") == 0) {

            if ( l_bcsr_sparse_block_size != 0) {

                if ( l_transpose == 1 ) {
                    matrix_t l_bcsr_input_matrix_transposed = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix_transposed, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                    matrix_t l_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);

                    l_rc = precond_cg(  l_precision,
                                        l_transpose,
                                        l_sparse_input_matrix_transposed->n, 
                                        X, 
                                        l_bcsr_input_matrix_transposed,
                                        l_iM_transposed,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                } else {
                    matrix_t l_bcsr_input_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                    matrix_t l_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);

                    l_rc = precond_cg(  l_precision,
                                        l_transpose,
                                        l_sparse_input_matrix->n, 
                                        X, 
                                        l_bcsr_input_matrix,
                                        l_iM,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                }
            } else {
                if ( l_transpose == 1 ) {
                    matrix_t l_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);

                    l_rc = precond_cg(  l_precision,
                                        l_transpose,
                                        l_sparse_input_matrix_transposed->n, 
                                        X, 
                                        l_sparse_input_matrix_transposed,
                                        l_iM_transposed,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                } else {
                    matrix_t l_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);

                    l_rc = precond_cg(  l_precision,
                                        l_transpose,
                                        l_sparse_input_matrix->n, 
                                        X, 
                                        l_sparse_input_matrix,
                                        l_iM,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                }
            }
        } else if (strcmp(l_solver_name, "QMR") == 0) {            
            if ( l_bcsr_sparse_block_size != 0) {
                matrix_t l_bcsr_input_matrix = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);
                matrix_t l_bcsr_input_matrix_transposed = VPFloatPackage::OSKIHelper::toBCSR(l_oski_sparse_input_matrix_transposed, l_bcsr_sparse_block_size, l_bcsr_sparse_block_size);

                l_rc = qmr(l_precision,
                            l_transpose,
                            l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                            X, 
                            l_transpose == 1 ? l_bcsr_input_matrix_transposed : l_bcsr_input_matrix, 
                            l_transpose == 1 ? l_bcsr_input_matrix : l_bcsr_input_matrix_transposed, 
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
            } else {
                l_rc = qmr(l_precision,
                            l_transpose,
                            l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                            X, 
                            l_transpose == 1 ? l_sparse_input_matrix_transposed : l_sparse_input_matrix, 
                            l_transpose == 1 ? l_sparse_input_matrix : l_sparse_input_matrix_transposed,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
            }
        } else {
            printf("Solver %s is not supported.\n", l_solver_name);
            exit(1);
        }

    } else {
        /*
         * DENSE version of solvers
         */
        matrix_t l_dense_input_matrix = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix, false, l_lda);

        if ( strcmp(l_solver_name, "BICG") == 0 || strcmp(l_solver_name, "BICGSTAB") == 0 || strcmp(l_solver_name, "PRECOND_BICG") == 0) {
            matrix_t l_dense_input_matrix_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix_transposed, false, l_lda);
            if ( strcmp(l_solver_name, "BICG") == 0 ) {
                l_rc = bicg(l_precision, 
                            l_transpose,
                            l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                            X, 
                            l_transpose == 1 ? l_dense_input_matrix_transposed : l_dense_input_matrix, 
                            l_transpose == 1 ? l_dense_input_matrix : l_dense_input_matrix_transposed, 
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
            } else if ( strcmp(l_solver_name, "BICGSTAB") == 0 ) {
                l_rc = bicgstab(l_precision, 
                                l_transpose,
                                l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                                X, 
                                l_transpose == 1 ? l_dense_input_matrix_transposed : l_dense_input_matrix, 
                                l_transpose == 1 ? l_dense_input_matrix : l_dense_input_matrix_transposed, 
                                B, 
                                l_tolerance, 
                                l_exponent_size, 
                                l_stride_size, 
                                l_log_buffer, 
                                l_log_buffer_size);
            } else {
                if ( l_transpose == 1 ) {
                    matrix_t l_sparse_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);
                    oski_matrix_wrapper_t l_oski_sparse_iM_transposed = VPFloatPackage::OSKIHelper::fromCSRMatrix(l_sparse_iM_transposed);
                    matrix_t l_dense_iM_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_iM_transposed, false, l_lda);

                    l_rc = precond_bicg(l_precision, 
                            l_transpose,
                            l_sparse_input_matrix_transposed->n, 
                            X, 
                            l_dense_input_matrix_transposed, 
                            l_dense_input_matrix, 
                            l_dense_iM_transposed,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
                } else {
                        matrix_t l_sparse_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);
                        oski_matrix_wrapper_t l_oski_sparse_iM = VPFloatPackage::OSKIHelper::fromCSRMatrix(l_sparse_iM);
                        matrix_t l_dense_iM = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_iM, false, l_lda);

                        l_rc = precond_bicg(l_precision, 
                            l_transpose,
                            l_sparse_input_matrix->n, 
                            X, 
                            l_dense_input_matrix, 
                            l_dense_input_matrix_transposed, 
                            l_dense_iM,
                            B, 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
                }
            }

        } else if (strcmp(l_solver_name, "CG") == 0) {
            if ( l_transpose == 1 ) { 
                matrix_t l_dense_input_matrix_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix_transposed, false, l_lda);
                l_rc = cg(  l_precision, 
                            l_transpose,
                            l_sparse_input_matrix_transposed->n, 
                            X, 
                            l_dense_input_matrix_transposed, 
                            B , 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
            } else {
                l_rc = cg(  l_precision, 
                            l_transpose,
                            l_sparse_input_matrix->n, 
                            X, 
                            l_dense_input_matrix, 
                            B , 
                            l_tolerance, 
                            l_exponent_size, 
                            l_stride_size, 
                            l_log_buffer, 
                            l_log_buffer_size);
            }

        } else if (strcmp(l_solver_name, "PRECOND_CG") == 0) {

                if ( l_transpose == 1 ) {
                    matrix_t l_dense_input_matrix_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix_transposed, false, l_lda);
                    matrix_t l_sparse_iM_transposed = jacobi(l_sparse_input_matrix_transposed, l_jacobi_shifter);
                    oski_matrix_wrapper_t l_oski_sparse_iM_transposed = VPFloatPackage::OSKIHelper::fromCSRMatrix(l_sparse_iM_transposed);
                    matrix_t l_dense_iM_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_iM_transposed, false, l_lda);

                    l_rc = precond_cg(  l_precision,
                                        l_transpose,
                                        l_sparse_input_matrix_transposed->n, 
                                        X, 
                                        l_dense_input_matrix_transposed,
                                        l_dense_iM_transposed,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                } else {
                    matrix_t l_sparse_iM = jacobi(l_sparse_input_matrix, l_jacobi_shifter);
                    oski_matrix_wrapper_t l_oski_sparse_iM = VPFloatPackage::OSKIHelper::fromCSRMatrix(l_sparse_iM);
                    matrix_t l_dense_iM = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_iM, false, l_lda);

                    l_rc = precond_cg(  l_precision, 
                                        l_transpose,
                                        l_sparse_input_matrix->n, 
                                        X, 
                                        l_dense_input_matrix,
                                        l_dense_iM,
                                        B, 
                                        l_tolerance, 
                                        l_exponent_size, 
                                        l_stride_size,
                                        l_log_buffer, 
                                        l_log_buffer_size);
                } 
        } else if ( strcmp(l_solver_name, "QMR") == 0 ) {
            matrix_t l_dense_input_matrix_transposed = VPFloatPackage::OSKIHelper::toDense(l_oski_sparse_input_matrix_transposed, false, l_lda);

            l_rc = qmr( l_precision,
                        l_transpose,
                        l_transpose == 1 ? l_sparse_input_matrix_transposed->n : l_sparse_input_matrix->n, 
                        X, 
                        l_transpose == 1 ? l_dense_input_matrix_transposed : l_dense_input_matrix, 
                        l_transpose == 1 ? l_dense_input_matrix : l_dense_input_matrix_transposed, 
                        B, 
                        l_tolerance, 
                        l_exponent_size, 
                        l_stride_size, 
                        l_log_buffer, 
                        l_log_buffer_size);
        } else {
            printf("Solver %s is not supported.\n", l_solver_name);
            exit(1);
        }

    }

    if (l_rc < 0 ) {
        printf("%s call failed %d.\n", l_solver_name, l_rc);
        l_rc = 1;
    } else {
        printf("%s call done.\n", l_solver_name);
        l_rc = 0;
    }

    if ( l_log_buffer_size > 0 ) {
        printf("%s Logs BEGIN ==========\n", l_solver_name);
        printf("%s", l_log_buffer);
        printf("%s Logs END ==========\n", l_solver_name);

        free(l_log_buffer);
    }
    
    free(X);
    free(B);
    free(l_solver_name);
    free(l_matrix_file_path);
    freeMatrix(l_sparse_input_matrix);
    return l_rc;
}
