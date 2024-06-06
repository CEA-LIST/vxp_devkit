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
#ifndef __VRP_IRQ_H_
#define __VRP_IRQ_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define VRP_MSI_BASE_ADDR               0x20100200000ULL
#define VRP_MSI_MASK_IT_0_MB_0_OFFSET 	0x800
#define VRP_MSI_MASK(it,mb) 			VRP_MSI_MASK_IT_0_MB_0_OFFSET + 0x80*it + 0x8*mb
#define VRP_MSI_STATUS_IT_0_MB_0_OFFSET 0x1000
#define VRP_MSI_STATUS(it,mb)			VRP_MSI_STATUS_IT_0_MB_0_OFFSET + 0x80*it + 0x8*mb
#define VRP_MSI_MAILBOX(mb) 	        0x10000*(mb+1)
#define VRP_MSI_MAILBOX_SUB_IT(mb)      VRP_MSI_MAILBOX(mb) + 0x8000

#define VRP_MSI_IT                      0
#define VRP_MSI_MB                      2
#define VRP_MSI_MB_SUB_IT               4 

/* Functions */
void vrp_msi_it_send(uint64_t msg);

#endif /* __VRP_IRQ_H_ */
