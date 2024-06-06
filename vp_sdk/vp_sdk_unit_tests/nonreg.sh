#!/bin/bash
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

export SILENT_MODE=1

if [ -z "${PKG_CONFIG_PATH}" ]
then
    echo "Set PKG_CONFIG_PATH before calling this script"
    exit 1
fi

################################################################################
# Build command for the VRP environment
################################################################################
function build_vrp() {

    TEST_DIRECTORY=${1}

    if [ -f ${TEST_DIRECTORY}/Makefile.vrp ]
    then

        RISCV_COMPILER_INSTALLATION_PATH=/home/560.1.361-EPI/xtools/riscv64-unknown-elf-newlib/bin
        CC=${RISCV_COMPILER_INSTALLATION_PATH}/riscv64-unknown-elf-gcc
        CXX=${RISCV_COMPILER_INSTALLATION_PATH}/riscv64-unknown-elf-g++
        RISCV_COMPILER_VRP_ISA_SPECIFIC_INSTALLATION_PATH=/home/560.1.361-EPI/vrp-binutils/Scientific-7.9/isa-vrp-v1.1.0/bin/
        OC=${RISCV_COMPILER_VRP_ISA_SPECIFIC_INSTALLATION_PATH}/riscv64-unknown-elf-objcopy
        OD=${RISCV_COMPILER_VRP_ISA_SPECIFIC_INSTALLATION_PATH}/riscv64-unknown-elf-objdump
        BIN2MEM=/home/560.1.361-EPI/SHARE/VRP_SDK/tools/bin2mem.py

        make -C ${TEST_DIRECTORY}  BSP=fpga_vcu128 CC=${CC} CXX=${CXX} OC=${OC} OD=${OD} BIN2MEM=${BIN2MEM} -f Makefile.vrp clean && \
        make -C ${TEST_DIRECTORY} -j 10 BSP=fpga_vcu128 CC=${CC} CXX=${CXX} OC=${OC} OD=${OD} BIN2MEM=${BIN2MEM} -f Makefile.vrp mem && \
        make -C ${TEST_DIRECTORY} -j 10 BSP=fpga_vcu128 CC=${CC} CXX=${CXX} OC=${OC} OD=${OD} BIN2MEM=${BIN2MEM} -f Makefile.vrp dump 
    fi

    return ${?}
}

################################################################################
# Run command for the VRP environment
################################################################################
function run_vrp() {

    TEST_DIRECTORY=${1}

    if [ -f ${TEST_DIRECTORY}/Makefile.vrp ]
    then
        sudo sh -c "../../vrp_run/vrp_run ${TEST_DIRECTORY} ${TEST_DIRECTORY}.x.bin ${REDIRECT_OUTPUT}"
    fi

    return ${?}
}

################################################################################
# Build command for the LINUX environment
################################################################################
function build_linux() {

    TEST_DIRECTORY=${1}
    
    if [ -f ${TEST_DIRECTORY}/Makefile.linux ]
    then
        make -C ${TEST_DIRECTORY} -f Makefile.linux clean && \
        make -C ${TEST_DIRECTORY} -f Makefile.linux 
    fi
    
    return ${?}
}

################################################################################
# Run command for the LINUX environment
################################################################################
function run_linux() {

    TEST_DIRECTORY=${1}

    if [ -f ${TEST_DIRECTORY}/Makefile.linux ]
    then
        ./${TEST_DIRECTORY}/${TEST_DIRECTORY} ${REDIRECT_OUTPUT}
    fi

    return ${?}
}

function contains() {
    [[ $1 =~ (^|[[:space:]])$2($|[[:space:]]) ]] && return 0 || return 1
}

TEST_ENVIRONMENT="all"
EXCLUDE_LIST="" 

# Process command line arguments
while [ ${#} -ne 0 ]
do
    case ${1} in 
        "--linux")
            TEST_ENVIRONMENT="linux"
            shift
            ;;
        "--vrp")
            TEST_ENVIRONMENT="vrp"
            shift
            ;;
        "--verbose")
            SILENT_MODE=0
            shift
            ;;
        "--exclude")
            shift
            EXCLUDE_LIST="test_${1} ${EXCLUDE_LIST}"
            shift
            ;;
        "--matrix_repo_path")
            shift
            export MATRIX_REPO_PATH="${1}"
            shift
            ;;
        *)
            break
            ;;
    esac
done

# Is there a specific test specified ?
if [ ${#} -ge 1 ] 
then
    # Yes ! use it
    TEST_LIST=test_${1}
else
    # No! Get all available tests
    TEST_LIST=`ls -d test_*`
fi

if [ ${SILENT_MODE} -eq 1 ]
then
    export MAKEFLAGS="--silent"
    export REDIRECT_OUTPUT="2>&1 >/dev/null"
fi

RC=0
READY_TO_RUN_TEST_LIST=""
FAIL_COMPILATION_TEST_LIST=""
COMPILATION_ERROR_COUNT=0

# Iterate over desired test and build them
for CURRENT_TEST in ${TEST_LIST}
do
    contains "${EXCLUDE_LIST}" "${CURRENT_TEST}"

    if [ ${?} -eq 0 ] 
    then
        echo "Skipping ${CURRENT_TEST} since it is excluded"
        continue
    fi

    if [ -d ${CURRENT_TEST} ]
    then
        case ${TEST_ENVIRONMENT} in
            "linux")
                build_linux ${CURRENT_TEST}
                ;;
            "vrp")
                build_vrp ${CURRENT_TEST}
                ;;
            "all")
                build_linux ${CURRENT_TEST} && \
                build_vrp ${CURRENT_TEST}
                ;;
        esac

        if [ ${?} -ne 0 ] 
        then
            COMPILATION_ERROR_COUNT=$((${COMPILATION_ERROR_COUNT} + 1))
            FAIL_COMPILATION_TEST_LIST="${FAIL_COMPILATION_TEST_LIST} ${CURRENT_TEST}"
        else
            READY_TO_RUN_TEST_LIST="${READY_TO_RUN_TEST_LIST} ${CURRENT_TEST}"
        fi
    else
        echo "test ${CURRENT_TEST} does not exists"
        RC=1
    fi
done

RUN_ERROR_COUNT=0
SUCCESSFUL_TEST_LIST=""

# Iterate over available tests and run them
for CURRENT_TEST in ${READY_TO_RUN_TEST_LIST}
do
    if [ -d ${CURRENT_TEST} ]
    then
        LINUX_RUN_RC=0
        VRP_RUN_RC=0
        case ${TEST_ENVIRONMENT} in
            "linux")
                build_linux ${CURRENT_TEST} && \
                run_linux ${CURRENT_TEST}
                LINUX_RUN_RC=${?}
                ;;
            "vrp")
                build_vrp ${CURRENT_TEST} && \
                run_vrp ${CURRENT_TEST}
                VRP_RUN_RC=${?}
                ;;
            "all")
                build_linux ${CURRENT_TEST} && \
                run_linux ${CURRENT_TEST} 
                LINUX_RUN_RC=${?}
                build_vrp ${CURRENT_TEST} && \
                run_vrp ${CURRENT_TEST}
                VRP_RUN_RC=${?}
                ;;
        esac

        if [ $((${VRP_RUN_RC} + ${LINUX_RUN_RC})) -ne 0 ] 
        then
            RUN_ERROR_COUNT=$((${RUN_ERROR_COUNT} + 1))
            FAIL_RUN_TEST_LIST="${FAIL_RUN_TEST_LIST} ${CURRENT_TEST}"
        else
            SUCCESSFUL_TEST_LIST="${SUCCESSFUL_TEST_LIST} ${CURRENT_TEST}"
        fi
    else
        echo "test ${CURRENT_TEST} does not exists"
        RC=1
    fi
done

# Exit when compilation error occurs
if [ ${COMPILATION_ERROR_COUNT} -ne 0 ]
then
    echo "Test failing in compilation(${COMPILATION_ERROR_COUNT}): ${FAIL_COMPILATION_TEST_LIST}"
    RC=1
fi

# Exit when run error occurs
if [ ${RUN_ERROR_COUNT} -ne 0 ]
then
    echo "Test run failed(${RUN_ERROR_COUNT}): ${FAIL_RUN_TEST_LIST}"
    RC=1
fi

if [ ! -z "${SUCCESSFUL_TEST_LIST}" ]
then
    echo "Test run succeed for: ${SUCCESSFUL_TEST_LIST}"
fi

exit ${RC}
