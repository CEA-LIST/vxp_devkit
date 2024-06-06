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

#include "VPSDK/VPComplex.hpp"

namespace VPFloatPackage {

    VPComplex::VPComplex( void * a_data, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, bool a_release_memory_on_destruction) :
        m_data(a_data)
    {
        this->m_environment.es = a_exponent_size;
        this->m_environment.bis = a_bis;
        this->m_environment.stride = a_stride;
        this->m_release_m_data_on_destruction = a_release_memory_on_destruction;
    }

    VPComplex::~VPComplex() {
        if ( this->m_release_m_data_on_destruction ) {
            free(this->m_data);
            this->m_data = NULL;
        }
    }

    VPComplex operator+ (const VPComplex& a_lhs, const VPComplex& a_rhs) {
        VPComplex l_result(a_lhs); l_result += a_rhs; return l_result;
    }

    VPComplex operator- (const VPComplex& a_lhs, const VPComplex& a_rhs) {
        VPComplex l_result(a_lhs); l_result -= a_rhs; return l_result;
    }

    VPComplex operator* (const VPComplex& a_lhs, const VPComplex& a_rhs) {
        VPComplex l_result(a_lhs); l_result *= a_rhs; return l_result;
    }

    VPComplex operator* (const VPComplex& a_lhs, const VPFloat& a_rhs) {
        VPComplex l_result(a_lhs); l_result.real(l_result.real()* a_rhs); l_result.imag(l_result.imag()* a_rhs); return l_result;
    }

    VPComplex operator/ (const VPComplex& a_lhs, const VPComplex& a_rhs) {
        VPComplex l_result(a_lhs); l_result /= a_rhs; return l_result;
    }

}