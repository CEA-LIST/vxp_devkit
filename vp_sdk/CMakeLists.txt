# Copyright 2023 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# 
# Authors       : Jerome Fereyre
# Creation Date : August, 2023
# Description   : 

cmake_minimum_required(VERSION 3.16)

string(TOLOWER ${VRP_PLATFORM} VRP_LOWER_PLATFORM)

# set the project name
project(vp_sdk_${VRP_PLATFORM})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(VBLAS_ENABLE_PERFMONITOR 1)
else()
    set(VBLAS_ENABLE_PERFMONITOR 0)
endif()

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 

pkg_check_modules(MATRIX_SDK_PKG REQUIRED IMPORTED_TARGET matrix_sdk_${VRP_LOWER_PLATFORM})
pkg_check_modules(MPFR_PKG REQUIRED IMPORTED_TARGET mpfr)
pkg_check_modules(MPC_PKG REQUIRED IMPORTED_TARGET mpc)
pkg_check_modules(GMP_PKG REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(VRP_SDK_PKG REQUIRED IMPORTED_TARGET vrp_sdk_${VRP_LOWER_PLATFORM})

# Creation of a static library (archive)
add_library(${PROJECT_NAME} STATIC)

FILE(GLOB VP_SDK_SOURCES   src/solvers/*.cpp
                            src/SparseMatrix/mmio.c
                            src/SparseMatrix/crs.c
                            src/SparseMatrix/csrIO.c)

set(VP_SDK_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/include)

list (APPEND VP_SDK_SOURCES src/VPSolvers/bicg/bicg_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSolvers/bicgstab/bicgstab_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSolvers/precond_bicg/precond_bicg_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSolvers/cg/cg_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSolvers/precond_cg/precond_cg_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSolvers/qmr/qmr_kernel.cpp)
list (APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplex_common.cpp)
list (APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplexArray_common.cpp)
list (APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASConfig.cpp)
list (APPEND VP_SDK_SOURCES src/VRPOffload/vrp_Matrix_serializer.cpp)
list (APPEND VP_SDK_SOURCES src/VRPOffload/vrp_Matrix_DENSE_serializer.cpp)
list (APPEND VP_SDK_SOURCES src/VRPOffload/vrp_Matrix_CSR_serializer.cpp)
list (APPEND VP_SDK_SOURCES src/VRPOffload/vrp_Matrix_BCSR_serializer.cpp)

list(APPEND pc_req_public "matrix_sdk_${VRP_PLATFORM}")
list(APPEND pc_req_public "vrp_sdk_${VRP_PLATFORM}")

string(COMPARE EQUAL ${VRP_LOWER_PLATFORM} linux_x86_64 _cmp)
if ( _cmp )
    ##################################################################
    ## LINUX x86_64
    ##################################################################
    
    list(APPEND VP_SDK_SOURCES src/VRPOffload/vrp_argument.cpp)
    list(APPEND VP_SDK_SOURCES src/VRPOffload/vrp_argument_array.cpp)
    list(APPEND VP_SDK_SOURCES src/VRPOffload/vrp_driver_interface.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/bicg/bicg_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/precond_bicg/precond_bicg_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/bicgstab/bicgstab_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/cg/cg_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/precond_cg/precond_cg_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/qmr/qmr_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPFloatpp/VPFloat_MPFR.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLAS_MPFR.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASComplex_MPFR.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASConfig_Linux.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VMath/VMath_MPFR.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplex_MPFR.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplexArray_MPFR.cpp)

    target_link_libraries(${PROJECT_NAME}
        PRIVATE
            -static
            ${MPFR_PKG_LIBRARIES}
            ${MPC_PKG_LIBRARIES}
            ${GMP_PKG_LIBRARIES}
            stdc++
            m
    )

    target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${MPFR_PKG_INCLUDE_DIRS}
            ${MPC_PKG_INCLUDE_DIRS}
            ${GMP_PKG_INCLUDE_DIRS}
            ${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include
            ${VRP_SDK_PKG_INCLUDE_DIRS}
            ${MATRIX_SDK_PKG_INCLUDE_DIRS}
    )

    target_compile_options(${PROJECT_NAME}
        PRIVATE
            -DUSE_MPFR
            -DVBLAS_ENABLE_PERFMONITOR=${VBLAS_ENABLE_PERFMONITOR}
    )   

    target_compile_options(${PROJECT_NAME}
        PUBLIC
            -I${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include
    )

    list(APPEND pc_req_public "mpfr")
    list(APPEND pc_req_public "mpc")
    list(APPEND pc_req_public "gmp")

else()
    ##################################################################
    ## RISCV VRP Target Configuration
    ##################################################################
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPFloatpp/VPFloat_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLAS_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASComplex_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASConfig_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VMath/VMath_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplex_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSDK/VPComplex/VPComplexArray_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/bicg/bicg_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/precond_bicg/precond_bicg_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/bicgstab/bicgstab_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/cg/cg_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/precond_cg/precond_cg_VRP.cpp)
    list(APPEND VP_SDK_SOURCES src/VPSolvers/qmr/qmr_VRP.cpp)

    pkg_check_modules(VRP_RISCV_BARE_PKG REQUIRED IMPORTED_TARGET vrp_riscv_bare_${BSP})
    
    #################################################################
    ## Debug mode addition configuration
    #################################################################
    if (${VBLAS_ENABLE_PERFMONITOR} EQUAL 1)
        list (APPEND VP_SDK_SOURCES src/VPSDK/VBLAS/VBLASPerfMonitor.cpp)
    endif()

    target_link_libraries(${PROJECT_NAME}
        PRIVATE
            ${VRP_RISCV_BARE_PKG_LIBRARIES}
            c
            c++
            gcc
    )
    
    target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${VRP_RISCV_BARE_PKG_INCLUDE_DIRS}
            ${VRP_SDK_PKG_INCLUDE_DIRS}
            ${MATRIX_SDK_PKG_INCLUDE_DIRS}
    )
    
    target_compile_options(${PROJECT_NAME}
        PRIVATE
            -DBSP_CONFIG_NCPUS=${BSP_CONFIG_NCPUS}
            -DCORE_REFCLK=${CORE_REFCLK}
            -DVBLAS_ENABLE_PERFMONITOR=${VBLAS_ENABLE_PERFMONITOR}
    )   

    list(APPEND pc_req_public "vrp_riscv_bare_${BSP}")

endif()




#################################################################
## Common target configuration
#################################################################
target_sources(${PROJECT_NAME}
    PRIVATE
        ${VP_SDK_SOURCES}
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${VP_SDK_INCLUDE_PATH}
)

get_target_property(VP_SDK_C_COMPILE_OPTIONS ${PROJECT_NAME} COMPILE_OPTIONS)
list(JOIN VP_SDK_C_COMPILE_OPTIONS " " VP_SDK_C_COMPILE_OPTIONS)

################################################################
## Generation of the pkg-config file to use this library
################################################################
list(JOIN pc_req_public " " pc_req_public)

configure_file(vp_sdk.pc.in ${PROJECT_NAME}.pc @ONLY)

################################################################
## Configuration of the library installation process
################################################################
install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc DESTINATION pkgconfig
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h*"
)
