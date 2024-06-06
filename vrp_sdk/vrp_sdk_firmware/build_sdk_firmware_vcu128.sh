#!/bin/csh
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

setenv VRP_RISCV_BARE_PATH `readlink -f ../../../vrp_linux_driver/vrp_riscv_bare`

setenv DISABLE_NC 1
echo ${PWD}
cd ${VRP_RISCV_BARE_PATH}
source setup.csh
cd -

set VRP_DRIVER_INCLUDE_PATH=`readlink -f ../../vrp_driver/include`

cd ${VRP_RISCV_BARE_PATH}/lib/vrp-sdk
make RISCVBARECEA=1 BSP=fpga_vcu128 BSP_CONFIG_NCPUS=1 UART_REFCLK=83000000 UART_BAUDRATE=38400 VRP_DRIVER_INCLUDE_PATH=${VRP_DRIVER_INCLUDE_PATH} clean
make RISCVBARECEA=1 BSP=fpga_vcu128 BSP_CONFIG_NCPUS=1 UART_REFCLK=83000000 UART_BAUDRATE=38400 VRP_DRIVER_INCLUDE_PATH=${VRP_DRIVER_INCLUDE_PATH}
cd -

make  clean

foreach VRP_SDK_FIRMWARE ( bicg cg precond_cg qmr)

    make BUILD_DIR=build_vrp_vcu128 BSP=fpga_vcu128 BSP_CONFIG_NCPUS=1 UART_REFCLK=83000000 UART_BAUDRATE=38400 VRP_DATA_ADDRESS=0x0000800200000000ULL VRP_SDK_FIRMWARE=${VRP_SDK_FIRMWARE} VRP_DRIVER_INCLUDE_PATH=${VRP_DRIVER_INCLUDE_PATH}
    make BUILD_DIR=build_vrp_vcu128 BSP=fpga_vcu128 BSP_CONFIG_NCPUS=1 UART_REFCLK=83000000 UART_BAUDRATE=38400 VRP_DATA_ADDRESS=0x0000800200000000ULL VRP_SDK_FIRMWARE=${VRP_SDK_FIRMWARE} VRP_DRIVER_INCLUDE_PATH=${VRP_DRIVER_INCLUDE_PATH} mem

    if ( ${?} != 0 ) then
        exit 1
    endif
    
    readlink -f wrappers/build_vrp/vrp_solver_${VRP_SDK_FIRMWARE}_vrp.x.mem 

end
