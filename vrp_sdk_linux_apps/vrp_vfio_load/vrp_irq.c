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
#include "vrp_irq.h"

void vrp_msi_it_send(uint64_t msg) {
    uint64_t addr;
    uint64_t data;

    // select sub IT 4
    uint64_t sub_it = (1 << VRP_MSI_MB_SUB_IT);

    // route sub-IT 4 of mailbox 2 to output IT 0
    addr = VRP_MSI_BASE_ADDR + VRP_MSI_MASK(VRP_MSI_IT, VRP_MSI_MB);
    data = sub_it;
    //printf("%s : 0x%llx = 0x%llx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;

    // trigger sub-IT 4 of mailbox 2
    addr = VRP_MSI_BASE_ADDR + VRP_MSI_MAILBOX(VRP_MSI_MB) + (sub_it*8);
    data = msg;
    //printf("%s : 0x%llx = 0x%llx\n", __func__, addr, data);
    *(uint64_t *)(addr) = data;
}
