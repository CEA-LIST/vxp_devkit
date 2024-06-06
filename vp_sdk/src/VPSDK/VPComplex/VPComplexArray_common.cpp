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

VPComplexArray::VPComplexArray(VPComplexArray & a_other) :
	m_nb_elements(a_other.m_nb_elements)
{	
	this->m_environment.es = a_other.m_environment.es;
	this->m_environment.bis = a_other.m_environment.bis;
	this->m_environment.stride = a_other.m_environment.stride;

	this->m_data = a_other.m_data;

	this->m_release_m_data_on_destruction = false;
}

VPComplexArray::VPComplexArray(VPComplexArray * a_other) :
	m_nb_elements(a_other->m_nb_elements)
{	
	this->m_environment.es = a_other->m_environment.es;
	this->m_environment.bis = a_other->m_environment.bis;
	this->m_environment.stride = a_other->m_environment.stride;

	this->m_data = a_other->m_data;

	this->m_release_m_data_on_destruction = false;
}

VPComplexArray::~VPComplexArray() {
	if ( this->m_release_m_data_on_destruction ) {
		free(this->m_data);
		this->m_data = NULL;
	}
}

void * VPComplexArray::operator+(int nb_elements) const {
	return (void *)((unsigned long long)(this->m_data) + (nb_elements * this->scalarSize()));
}

VPComplex VPComplexArray::operator[](int nb_elements) const {
	return VPComplex((void *)((unsigned long long)(this->m_data) + (nb_elements * this->scalarSize())), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
}