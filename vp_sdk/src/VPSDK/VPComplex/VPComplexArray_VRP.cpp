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
#include <cstring>

using namespace VPFloatPackage;

VPComplexArray::VPComplexArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride,  int a_nb_elements) :
    m_nb_elements(a_nb_elements)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;

	uint64_t l_scalar_size = this->scalarSize();

    this->m_data = malloc(l_scalar_size * 2 * this->m_nb_elements);

	memset(this->m_data, 0, l_scalar_size * 2 * this->m_nb_elements);

	this->m_release_m_data_on_destruction = true;
}

VPComplexArray::VPComplexArray(double * a_other, int a_nb_elements):
	m_nb_elements(a_nb_elements)
{
	this->m_environment.es = VPFLOAT_EVP_DOUBLE.es;
	this->m_environment.bis = VPFLOAT_EVP_DOUBLE.bis;
	this->m_environment.stride = VPFLOAT_EVP_DOUBLE.stride;

	this->m_data = a_other;

	this->m_release_m_data_on_destruction = false;
}

VPComplexArray::VPComplexArray(float * a_other, int a_nb_elements):
	m_nb_elements(a_nb_elements)
{
	this->m_environment.es = VPFLOAT_EVP_FLOAT.es;
	this->m_environment.bis = VPFLOAT_EVP_FLOAT.bis;
	this->m_environment.stride = VPFLOAT_EVP_FLOAT.stride;

	this->m_data = a_other;

	this->m_release_m_data_on_destruction = false;
}

uint64_t VPComplexArray::scalarSize() const {
    return VPFLOAT_SIZEOF(this->m_environment);
}
