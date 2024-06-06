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
project(vrp_sdk_${VRP_PLATFORM})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(VRP_SDK_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/include)

set(VRP_VBLAS_ENABLE_HWPF 1)
set(VRP_SPVBLAS_ENABLE_HWPF 1)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list (APPEND VRP_SDK_SOURCES src/VRPSDK/vblas_perfmonitor.c)
    set(VBLAS_ENABLE_PERFMONITOR 1)
else()
    set(VBLAS_ENABLE_PERFMONITOR 0)
endif()

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 

# these calls create special `PkgConfig::<MODULE>` variables
pkg_check_modules(MATRIX_SDK_PKG REQUIRED IMPORTED_TARGET matrix_sdk_${VRP_LOWER_PLATFORM})


# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 

# these calls create special `PkgConfig::<MODULE>` variables
pkg_check_modules(MATRIX_SDK_PKG REQUIRED IMPORTED_TARGET matrix_sdk_${VRP_LOWER_PLATFORM})

string(COMPARE NOTEQUAL ${VRP_LOWER_PLATFORM} linux_x86_64 _cmp)
if ( _cmp )

    pkg_check_modules(VRP_RISCV_BARE_PKG REQUIRED IMPORTED_TARGET vrp_riscv_bare_${VRP_LOWER_PLATFORM})

    ##################################################################
    ## RISCV VRP Target Configuration
    ##################################################################

    # Creation of a static library (archive)
    add_library(${PROJECT_NAME} STATIC)

    file(GLOB VRP_PERFCOUNTERS_SOURCES src/VRPSDK/perfcounters/*.c)
    file(GLOB VRP_VBLAS_SOURCES src/VRPSDK/vblas/*.c)
    file(GLOB VRP_VBLAS_MT_SOURCES src/VRPSDK/vblas_mt/*.c)
    file(GLOB VRP_VUTILS_SOURCES src/VRPSDK/vutils/*.c)
    file(GLOB VRP_VMATH_SOURCES src/VRPSDK/vmath/*.c)
    file(GLOB_RECURSE VRP_SPVBLAS_SOURCES src/VRPSDK/spvblas/*.c)
    list(APPEND VRP_SDK_SOURCES ${VRP_PERFCOUNTERS_SOURCES})
    list(APPEND VRP_SDK_SOURCES ${VRP_VBLAS_SOURCES})
    list(APPEND VRP_SDK_SOURCES ${VRP_VBLAS_MT_SOURCES})
    list(APPEND VRP_SDK_SOURCES ${VRP_VUTILS_SOURCES})
    list(APPEND VRP_SDK_SOURCES ${VRP_VMATH_SOURCES})
    list(APPEND VRP_SDK_SOURCES ${VRP_SPVBLAS_SOURCES})


    #################################################################
    ## Debug mode addition configuration
    #################################################################
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        list (APPEND VRP_SDK_SOURCES src/VRPSDK/vblas_perfmonitor.c)
        target_compile_options(${PROJECT_NAME}
            PRIVATE
                -DVBLAS_ENABLE_PERFMONITOR=${VBLAS_ENABLE_PERFMONITOR}
        )
    endif()

    target_link_libraries(${PROJECT_NAME}
        PRIVATE
            vrp_riscv_bare_${BSP}
            c
            c++
            gcc
    )
    
    target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${MATRIX_SDK_PKG_INCLUDE_DIRS}
            ${VRP_RISCV_BARE_PKG_INCLUDE_DIRS}
    )

    #################################################################
    ## Common target configuration
    #################################################################
    target_sources(${PROJECT_NAME}
        PRIVATE
            ${VRP_SDK_SOURCES}
    )

    target_include_directories(${PROJECT_NAME}
        PUBLIC
            ${VRP_SDK_INCLUDE_PATH}
            ${MATRIX_SDK_PKG_INCLUDE_DIRS}
    )

    target_compile_options(${PROJECT_NAME}
        PUBLIC
            -DVBLAS_ENABLE_HWPF=${VRP_VBLAS_ENABLE_HWPF}
            -DSPVBLAS_ENABLE_HWPF=${VRP_SPVBLAS_ENABLE_HWPF}
            -DVBLAS_ENABLE_PERFMONITOR=${VBLAS_ENABLE_PERFMONITOR}
    )

    ################################################################
    ## Configuration of the library installation process
    ################################################################
    install(TARGETS ${PROJECT_NAME}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib/static
    )

    list(APPEND pc_req_public "vrp_riscv_bare_${VRP_PLATFORM}")

    get_target_property(VRP_SDK_C_COMPILE_OPTIONS ${PROJECT_NAME} COMPILE_OPTIONS)
    list(JOIN VRP_SDK_C_COMPILE_OPTIONS " " VRP_SDK_C_COMPILE_OPTIONS)

    ################################################################
    ## Generation of the pkg-config file to use this library
    ################################################################
    list(APPEND pc_req_public "matrix_sdk_${VRP_PLATFORM}")

    list(JOIN pc_req_public " " pc_req_public)

    configure_file(vrp_sdk.pc.in ${PROJECT_NAME}.pc @ONLY)

else()

    ################################################################
    ## Generation of the pkg-config file to use this library
    ################################################################
    list(APPEND pc_req_public "matrix_sdk_${VRP_PLATFORM}")

    list(JOIN pc_req_public " " pc_req_public)

    configure_file(vrp_sdk_headers_only.pc.in ${PROJECT_NAME}.pc @ONLY)

endif()

################################################################
## Configuration of the library installation process
################################################################
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc DESTINATION pkgconfig
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h*"
)