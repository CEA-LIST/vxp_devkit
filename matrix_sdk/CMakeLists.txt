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
# Authors       : Jerome Fereyre
# Creation Date : August, 2023
# Description   : 

cmake_minimum_required(VERSION 3.16)

string(TOLOWER ${VRP_PLATFORM} VRP_LOWER_PLATFORM)

# set the project name
project(matrix_sdk_${VRP_PLATFORM})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 

# Creation of a static library (archive)
add_library(${PROJECT_NAME} STATIC)

set(MATRIX_SDK_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/include)

list (APPEND MATRIX_SDK_SOURCES src/MTXUtil/mmio.c)
list (APPEND MATRIX_SDK_SOURCES src/MTXUtil/crs.c)
list (APPEND MATRIX_SDK_SOURCES src/MTXUtil/crsIO.c)
list (APPEND MATRIX_SDK_SOURCES src/Matrix/matrix.cpp)
list (APPEND MATRIX_SDK_SOURCES src/Matrix/matrixComplex.cpp)

string(COMPARE EQUAL ${VRP_LOWER_PLATFORM} linux_x86_64 _cmp)
if ( _cmp )
    ##################################################################
    ## LINUX x86_64
    ##################################################################
    
    list(APPEND MATRIX_SDK_SOURCES src/MTXUtil/MTXParser.cpp)
    list(APPEND MATRIX_SDK_SOURCES src/OSKIHelper/OSKIHelper.cpp)
    list(APPEND MATRIX_SDK_SOURCES src/Preconditionners/jacobi.cpp)

    list (APPEND MATRIX_SDK_HEADERS ${PROJECT_SOURCE_DIR}/include/MTXUtil/MTXParser.hpp)
    list (APPEND MATRIX_SDK_HEADERS ${PROJECT_SOURCE_DIR}/include/OSKIHelper.hpp)
    list (APPEND MATRIX_SDK_HEADERS ${PROJECT_SOURCE_DIR}/include/Preconditionners.hpp)

    pkg_check_modules(OSKI_PKG REQUIRED IMPORTED_TARGET oski)

    target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${OSKI_PKG_INCLUDE_DIRS}
    )

    list(APPEND pc_req_public "oski")
endif()

#################################################################
## Common target configuration
#################################################################
target_sources(${PROJECT_NAME}
    PRIVATE
        ${MATRIX_SDK_SOURCES}
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${MATRIX_SDK_INCLUDE_PATH}
)

################################################################
## Generation of the pkg-config file to use this library
################################################################
configure_file(matrix_sdk.pc.in ${PROJECT_NAME}.pc @ONLY)

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
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/Matrix DESTINATION include FILES_MATCHING PATTERN "*.h*"
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/MTXUtil DESTINATION include FILES_MATCHING PATTERN "*.h*"
)

install(
    FILES ${MATRIX_SDK_HEADERS} DESTINATION include
)
