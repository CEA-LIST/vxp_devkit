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

#include <mpc.h>

using namespace VPFloatPackage;

VPComplexArray::VPComplexArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride,  int a_nb_elements) :
    m_nb_elements(a_nb_elements)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;

	uint64_t l_scalar_size = this->scalarSize();

    this->m_data = malloc(l_scalar_size * this->m_nb_elements);

	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	for (int l_index = 0 ; l_index < a_nb_elements; l_index++) {
		mpc_init2(((mpc_t *)this->m_data)[l_index], this->m_environment.bis - this->m_environment.es - 1 + 1);
		mpc_set_d_d(((mpc_t *)this->m_data)[l_index], 0.0, 0.0, mpfr_get_default_rounding_mode());
	}

	this->m_release_m_data_on_destruction = true;
}

VPComplexArray::VPComplexArray(double * a_other, int a_nb_elements):
	m_nb_elements(a_nb_elements)
{
	this->m_environment.es = VPFLOAT_EVP_DOUBLE.es;
	this->m_environment.bis = VPFLOAT_EVP_DOUBLE.bis;
	this->m_environment.stride = VPFLOAT_EVP_DOUBLE.stride;

	uint64_t l_scalar_size = this->scalarSize();

	this->m_data = malloc(l_scalar_size * this->m_nb_elements);

	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	for (int l_index = 0 ; l_index < a_nb_elements; l_index++) {

		mpc_init2(((mpc_t *)this->m_data)[l_index], this->m_environment.bis - this->m_environment.es - 1 + 1);
		mpc_set_d_d(((mpc_t *)this->m_data)[l_index], a_other[l_index*2], a_other[l_index*2+1], mpfr_get_default_rounding_mode());
	}

	this->m_release_m_data_on_destruction = true;
}

uint64_t VPComplexArray::scalarSize() const {
    return sizeof(mpc_t);
}