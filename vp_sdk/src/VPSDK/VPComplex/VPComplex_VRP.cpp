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
#include "VPSDK/VMath.hpp"
#include "VRPSDK/vmath.h"
#include <cstring>

using namespace VPFloatPackage;

void VPComplex::init_from_environment() {
	uint64_t l_scalar_size = VPFLOAT_SIZEOF(this->m_environment);

	this->m_data = malloc(l_scalar_size * 2);

	memset(this->m_data, 0, l_scalar_size * 2);
	
	this->m_release_m_data_on_destruction = true;
}

VPComplex::VPComplex() {
	this->m_environment.es = VPFloatComputingEnvironment::get_temporary_var_environment().es;
	this->m_environment.bis = VPFloatComputingEnvironment::get_temporary_var_environment().bis;
	this->m_environment.stride =  VPFloatComputingEnvironment::get_temporary_var_environment().stride;

	init_from_environment();
}

/*
 *
 */
VPComplex::VPComplex(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride) {
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;

	init_from_environment();
}

/*
 *
 */
VPComplex::VPComplex(const VPComplex &  a_other) {
	this->m_environment.es = a_other.m_environment.es;
	this->m_environment.bis = a_other.m_environment.bis;
	this->m_environment.stride = a_other.m_environment.stride;

	init_from_environment();

	VPFloat l_real(this->m_data, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
	VPFloat l_imag((void *)(((uint64_t)(this->m_data)) + VPFLOAT_SIZEOF(this->m_environment)), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	l_real = a_other.real();
	l_imag = a_other.imag();
}

const VPFloat VPComplex::real() const {
	return VPFloat(this->m_data, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
}

const VPFloat VPComplex::imag() const {
	return VPFloat((void *)(((uint64_t)(this->m_data)) + VPFLOAT_SIZEOF(this->m_environment)), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
}

void VPComplex::real(double a_real) {
	VPFloat l_real(this->m_data, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	l_real = a_real;
}

void VPComplex::imag(double a_imag) {
	VPFloat l_imag((void *)(((uint64_t)(this->m_data)) + VPFLOAT_SIZEOF(this->m_environment)), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	l_imag = a_imag;
}

VPComplex VPComplex::conjugate() const {
	VPComplex l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	::vassign(this->real().getData(), this->real().getEnvironment(), l_result.real().getData(), l_result.real().getEnvironment());

	::vassign(this->imag().getData(), this->imag().getEnvironment(), l_result.imag().getData(), l_result.imag().getEnvironment());
	::vinvertsign(VPFloatComputingEnvironment::get_precision(), l_result.imag().getData(), l_result.imag().getEnvironment());

	return l_result;
}

VPFloat VPComplex::norm() const {
	VPFloat l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	::vnorm(VPFloatComputingEnvironment::get_precision(), this->m_data, this->m_environment, l_result.getData(), l_result.getEnvironment());

	return l_result;
}

void VPComplex::add(const VPComplex & a_other) {

}

void VPComplex::sub(const VPComplex & a_other) {

}

void VPComplex::mul(const VPComplex & a_other) {

}

void VPComplex::div(const VPComplex & a_other) {

}

VPComplex& VPComplex::operator  = (const int x) {
	VPFloat l_real(this->m_data, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	l_real = x;

	return *this;
}

VPComplex& VPComplex::operator  = (const VPComplex &a_other) {
	VPFloat l_real(this->m_data, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
	VPFloat l_imag((void *)(((uint64_t)(this->m_data)) + VPFLOAT_SIZEOF(this->m_environment)), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	l_real = a_other.real();
	l_imag = a_other.imag();
	
	return *this;
}

namespace VPFloatPackage {

	bool operator==( const VPComplex & a_lhs, const VPComplex & a_rhs ) {
		return a_lhs.real() == a_rhs.real() && a_lhs.imag() == a_rhs.imag();
	}

	bool operator!=( const VPComplex & a_lhs, const VPComplex & a_rhs ) {
		return ! (a_lhs == a_rhs);
	}

	std::ostream& operator << (std::ostream & stream, const VPComplex & a_complex) {
		stream << "real:" << a_complex.real() << " - imag:" << a_complex.imag() ;
		return stream;
	}

}