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
#include <stdio.h>
#include <mpc.h>
#include <mpfr.h>
#include <cstring>

using namespace VPFloatPackage;

void VPComplex::init_from_environment() {
    this->m_data = malloc(sizeof(mpc_t));

	memset(this->m_data, 0, sizeof(mpc_t));

	this->m_release_m_data_on_destruction = true;
	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	mpc_init2(*(mpc_t *)(this->m_data),  this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpc_set_d_d(*(mpc_t *)(this->m_data), 0.0, 0.0, mpfr_get_default_rounding_mode());
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

	mpc_set(*(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data), mpfr_get_default_rounding_mode());
}

const VPFloat VPComplex::real() const {
	VPFloat l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
	
	mpc_real(*(mpfr_t *)(l_result.getData()), *(mpc_t *)(this->m_data), mpfr_get_default_rounding_mode());

	return l_result;
}

const VPFloat VPComplex::imag() const {
	VPFloat l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);

	mpc_imag(*(mpfr_t *)(l_result.getData()), *(mpc_t *)(this->m_data), mpfr_get_default_rounding_mode());

	return l_result;
}

void VPComplex::real(double a_real) {
	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	mpfr_set_d(mpc_realref((*(mpc_t *)(this->m_data))), a_real, mpfr_get_default_rounding_mode());
}

void VPComplex::imag(double a_imag) {
	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	mpfr_set_d(mpc_imagref(*(mpc_t *)(this->m_data)), a_imag, mpfr_get_default_rounding_mode());
}

VPComplex VPComplex::conjugate() const {
	VPComplex l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
	
	mpc_conj(*(mpc_t *)(l_result.m_data), *(mpc_t *)(this->m_data), mpfr_get_default_rounding_mode());
	return l_result;
}

VPFloat VPComplex::norm() const {
	VPFloat l_result(this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
	
	mpc_norm(*(mpfr_t *)(l_result.getData()), *(mpc_t *)(this->m_data), mpfr_get_default_rounding_mode());

	return l_result;
}

void VPComplex::add(const VPComplex & a_other) {
	mpc_add(*(mpc_t *)(this->m_data), *(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data), mpfr_get_default_rounding_mode());
}

void VPComplex::sub(const VPComplex & a_other) {
	mpc_sub(*(mpc_t *)(this->m_data), *(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data), mpfr_get_default_rounding_mode());
}

void VPComplex::mul(const VPComplex & a_other) {
	mpc_mul(*(mpc_t *)(this->m_data), *(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data), mpfr_get_default_rounding_mode());
}

void VPComplex::div(const VPComplex & a_other) {
	mpc_div(*(mpc_t *)(this->m_data), *(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data), mpfr_get_default_rounding_mode());
}

VPComplex& VPComplex::operator  = (const int x) {
	
	mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

	mpc_set_si(*(mpc_t *)(this->m_data), (long int)x, mpfr_get_default_rounding_mode());
	return *this;
}

VPComplex& VPComplex::operator  = (const VPComplex &a_other) {
	mpc_set(*(mpc_t *)(this->m_data), *(mpc_t *)(a_other.m_data) , mpfr_get_default_rounding_mode());
	return *this;
}

namespace VPFloatPackage {

	bool operator==( const VPComplex & a_lhs, const VPComplex & a_rhs ) {
		return (mpc_cmp(*(mpc_t *)(a_lhs.m_data), *(mpc_t *)(a_rhs.m_data)) == 0);
	}

	bool operator!=( const VPComplex & a_lhs, const VPComplex & a_rhs ) {
		return ! (a_lhs == a_rhs);
	}

	std::ostream& operator << (std::ostream & stream, const VPComplex & a_complex) {
		char * l_buffer = mpc_get_str(10, 0, *((mpc_t *)(a_complex.m_data)), mpfr_get_default_rounding_mode());
		
		stream << l_buffer << std::flush;

		mpc_free_str(l_buffer);
			
		return stream;


	}

}