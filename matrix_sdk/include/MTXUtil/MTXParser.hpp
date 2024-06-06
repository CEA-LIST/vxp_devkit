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

#ifndef __MTXPARSER_HPP__
#define __MTXPARSER_HPP__

#include "MTXUtil/mmio.h"
#include "MTXUtil/crs.h"
#include <string>

namespace MTXParser {

    /*
     * Function used to parse and MTX file and generate a CSR matrix.
     * Real an Complex values are supported.
     * Indeces in matrix are 1-based
     */
    crs_t* parseFileToCRS(std::string a_mtx_file_path);

    /*
     * Display CSR matrice with 1-based indeces.
     */
    void displayCSR(crs_t * SPM);

};

#endif /* __MTXPARSER_HPP__ */