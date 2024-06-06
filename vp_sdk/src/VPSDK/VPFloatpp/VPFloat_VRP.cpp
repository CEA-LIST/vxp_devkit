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

#include "VPSDK/VPFloat.hpp"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstddef>
#include "VRPSDK/vmath.h"
#include "VRPSDK/vutils.h"

using namespace VPFloatPackage;

/*******************************************************************************************************************
 * VPFloatComputingEnvironment
 ******************************************************************************************************************/
void VPFloatComputingEnvironment::set_rounding_mode(VPFloatRoundingMode a_rounding_mode) {
	vpfloat_set_rm_comp(a_rounding_mode);
}

uint8_t VPFloatComputingEnvironment::get_rounding_mode() {
	return vpfloat_get_rm_comp();
}


/*******************************************************************************************************************
 * VPFloat Class
 ******************************************************************************************************************/
VPFloat::VPFloat() :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = VPFloatComputingEnvironment::get_temporary_var_environment().es;
	this->m_environment.bis = VPFloatComputingEnvironment::get_temporary_var_environment().bis;
	this->m_environment.stride =  VPFloatComputingEnvironment::get_temporary_var_environment().stride;

	this->m_data = malloc(VPFLOAT_SIZEOF(this->m_environment));
	this->m_release_m_data_on_destruction = true;

	memset(this->m_data, 0, VPFLOAT_SIZEOF(this->m_environment));
}

VPFloat::VPFloat(const VPFloat & a_other) :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_other.m_environment.es;
	this->m_environment.bis = a_other.m_environment.bis;
	this->m_environment.stride = a_other.m_environment.stride;

	this->m_data = malloc(VPFLOAT_SIZEOF(this->m_environment));
	this->m_release_m_data_on_destruction = true;

	memcpy(this->m_data, a_other.m_data, VPFLOAT_SIZEOF(this->m_environment));
}

VPFloat::VPFloat(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride) :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;
	
	this->m_data = malloc(VPFLOAT_SIZEOF(this->m_environment));
	this->m_release_m_data_on_destruction = true;

	memset(this->m_data, 0, VPFLOAT_SIZEOF(this->m_environment));
}

VPFloat::VPFloat( void * a_data, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, bool a_release_memory_on_destruction) :
	m_data(a_data),
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;
	this->m_release_m_data_on_destruction = a_release_memory_on_destruction;
}

VPFloat::VPFloat( uint64_t a_mantissa, uint64_t a_exponent, bool a_sign, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride ) :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;

	this->m_data = NULL;

	this->m_data = ::vbuild(a_mantissa, a_exponent, a_sign, NULL, this->m_environment);
	this->m_release_m_data_on_destruction = true;
}

VPFloat::VPFloat( double a_value ) :
	m_double_data(true),
	m_float_data(false)
{
	this->m_environment.es = VPFLOAT_EVP_DOUBLE.es;
	this->m_environment.bis = VPFLOAT_EVP_DOUBLE.bis;
	this->m_environment.stride = VPFLOAT_EVP_DOUBLE.stride;

	this->m_data = malloc(VPFLOAT_SIZEOF(VPFLOAT_EVP_DOUBLE));
	this->m_release_m_data_on_destruction = true;

	::dtov(a_value, this->m_data, this->m_environment);
}

VPFloat::VPFloat( float a_value ) :
	m_double_data(false),
	m_float_data(true)
{
	this->m_environment.es = VPFLOAT_EVP_FLOAT.es;
	this->m_environment.bis = VPFLOAT_EVP_FLOAT.bis;
	this->m_environment.stride = VPFLOAT_EVP_FLOAT.stride;

	this->m_data = malloc(VPFLOAT_SIZEOF(VPFLOAT_EVP_FLOAT));
	this->m_release_m_data_on_destruction = true;

	::ftov(a_value, this->m_data, this->m_environment);
}

VPFloat::~VPFloat() {
	if ( this->m_release_m_data_on_destruction && this->m_data != NULL ) {
		free(this->m_data);
	}
}

vpfloat_evp_t VPFloat::getEnvironment() const {
	return this->m_environment;
}

void * VPFloat::getData() const {
	return this->m_data;
}

double VPFloat::getdouble() const {
	return ::vtod(this->m_data, this->m_environment);
}

float VPFloat::getfloat() const {
	return ::vtof(this->m_data, this->m_environment);
}

int64_t VPFloat::exponent() const {
	return ::vexponent(this->m_data, this->m_environment);
}

void VPFloat::exponent(const int64_t a_exponent) {
	::vexponent_set(this->m_data, this->m_environment, (uint64_t)a_exponent);
}

uint64_t VPFloat::mantissaChunk(const uint16_t a_chunk_index) const {
    switch (a_chunk_index) {
        case 0 :
            return ::vmantissaChunk(0, this->m_data, this->m_environment);
            break;
        case 1 :
            return ::vmantissaChunk(1, this->m_data, this->m_environment);
            break;
        case 2 :
            return ::vmantissaChunk(2, this->m_data, this->m_environment);
            break;
        case 3 :
            return ::vmantissaChunk(3, this->m_data, this->m_environment);
            break;
        case 4 :
            return ::vmantissaChunk(4, this->m_data, this->m_environment);
            break;
        case 5 :
            return ::vmantissaChunk(5, this->m_data, this->m_environment);
            break;
        case 6 :
            return ::vmantissaChunk(6, this->m_data, this->m_environment);
            break;
        case 7 :
            return ::vmantissaChunk(7, this->m_data, this->m_environment);
            break;
        default :
            std::cout << __func__ << " Invalid Chunk Index." << std::endl;
            return 0;
    }
}

 void VPFloat::mantissaChunk(const uint16_t a_chunk_index, const uint64_t a_chunk_value) {
	switch (a_chunk_index) {
        case 0 :
            ::vmantissaChunk_set(0, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 1 :
            ::vmantissaChunk_set(1, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 2 :
            ::vmantissaChunk_set(2, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 3 :
            ::vmantissaChunk_set(3, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 4 :
            ::vmantissaChunk_set(4, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 5 :
            ::vmantissaChunk_set(5, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 6 :
            ::vmantissaChunk_set(6, this->m_data, this->m_environment, a_chunk_value);
            break;
        case 7 :
            ::vmantissaChunk_set(7, this->m_data, this->m_environment, a_chunk_value);
            break;
        default :
            std::cout << __func__ << " Invalid Chunk Index." << std::endl;
            break;
    }
}

bool VPFloat::isNaN() const {
	return ( ::isNaN(this->m_data, this->m_environment) == 1 );
}

bool VPFloat::isInf() const {
	return ( ::isInf(this->m_data, this->m_environment) == 1 );
}

namespace VPFloatPackage {

	VPFloat abs(const VPFloat & a_x) {
		void * l_result = ::vabs(	VPFloatComputingEnvironment::get_precision(),
								a_x.m_data, a_x.m_environment,
								NULL, a_x.m_environment);
		return VPFloat (l_result, a_x.m_environment.es, a_x.m_environment.bis, a_x.m_environment.stride, true );
	}

	VPFloat operator+(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vadd( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								a_rhs.m_data, a_rhs.m_environment,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );
	}

	VPFloat operator+( const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vadd( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								(void *)&a_rhs, VPFLOAT_EVP_DOUBLE,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );
	}

	VPFloat operator-(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result_data = ::vsub( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								a_rhs.m_data, a_rhs.m_environment,
								NULL, l_tmp_var_environment);
		VPFloat l_result(l_result_data, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );
		return l_result;
	}

	VPFloat operator-(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vsub( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								(void *)&a_rhs, VPFLOAT_EVP_DOUBLE,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );	
	}

	VPFloat operator-(double a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vsub( VPFloatComputingEnvironment::get_precision(),
								(void *)&a_lhs, VPFLOAT_EVP_DOUBLE,
								a_rhs.m_data, a_rhs.m_environment,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );	
	}

	VPFloat operator*(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vmul( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								a_rhs.m_data, a_rhs.m_environment,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );
	}

	VPFloat operator*(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vmul( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								(void *)&a_rhs, VPFLOAT_EVP_DOUBLE,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );	
	}

	VPFloat operator/(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		// printf("== div == \n");
		// printf("div param a_lhs es=%d - bis=%d - stride=%d: double_conversion=%e - vpfloat_dump=%s\n", a_lhs.m_environment.es, a_lhs.m_environment.bis, a_lhs.m_environment.stride, ::vtod(a_lhs.m_data, a_lhs.m_environment), ::vtostring(a_lhs.m_data, a_lhs.m_environment));
		// printf("div param a_rhs es=%d - bis=%d - stride=%d: double_conversion=%e - vpfloat_dump=%s\n", a_rhs.m_environment.es, a_rhs.m_environment.bis, a_rhs.m_environment.stride, ::vtod(a_rhs.m_data, a_rhs.m_environment), ::vtostring(a_rhs.m_data, a_rhs.m_environment));
		// printf("div desired result environment es=%d - bis=%d - stride=%d.\n", l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		void * l_result = ::vdiv(	VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								a_rhs.m_data, a_rhs.m_environment,
								NULL, l_tmp_var_environment);
		return VPFloat (l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );
		
		
		// printf("1- div result with env es=%d - bis=%d - stride=%d: double_conversion=%e - vpfloat_dump=%s\n", l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, ::vtod(l_result, l_tmp_var_environment), ::vtostring(l_result, l_tmp_var_environment));
		
		// l_tmp_var_environment.es = 8;

		// printf("div desired result environment es=%d - bis=%d - stride=%d.\n", l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		// void * l_result2 = ::vdiv(	VPFloatComputingEnvironment::get_precision(),
		// 						a_lhs.m_data, a_lhs.m_environment,
		// 						a_rhs.m_data, a_rhs.m_environment,
		// 						NULL, l_tmp_var_environment);
		// VPFloat l_tmp_value2(l_result2, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );

		// printf("2- div result with env es=%d - bis=%d - stride=%d: double_conversion=%e - vpfloat_dump=%s\n", l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, ::vtod(l_result2, l_tmp_var_environment), ::vtostring(l_result2, l_tmp_var_environment));

		// return l_tmp_value;
	}

	VPFloat operator/(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		void * l_result = ::vdiv( VPFloatComputingEnvironment::get_precision(),
								a_lhs.m_data, a_lhs.m_environment,
								(void *)&a_rhs, VPFLOAT_EVP_DOUBLE,
								NULL, l_tmp_var_environment);
		return VPFloat(l_result, l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride, true );		
	}

	VPFloat operator-(const VPFloat & a_rhs) {
		void * l_result = ::vinvertsign(VPFloatComputingEnvironment::get_precision(), a_rhs.m_data, a_rhs.m_environment);
		return VPFloat(l_result, a_rhs.m_environment.es, a_rhs.m_environment.bis, a_rhs.m_environment.stride, true );				
	}

	bool operator>( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_gt(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	bool operator>=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_geq(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	bool operator<( const VPFloat & a_lhs, double a_rhs ) {
		VPFloat l_rhs(a_rhs);
		uint64_t l_result = ::vcompare_lt(a_lhs.m_data, a_lhs.m_environment, l_rhs.m_data, l_rhs.m_environment);
		return l_result;
	}

	bool operator<( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_lt(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	bool operator<=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_leq(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	bool operator==( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_eq(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	bool operator!=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return ::vcompare_neq(a_lhs.m_data, a_lhs.m_environment, a_rhs.m_data, a_rhs.m_environment);
	}

	std::ostream& operator << (std::ostream & stream, const VPFloat & a_vpfloat) {
		char * l_buffer = ::vtostring(a_vpfloat.m_data, a_vpfloat.m_environment);
		stream << l_buffer << std::flush;
		free(l_buffer);
		return stream;
	}
}

VPFloat & VPFloat::operator+=(const VPFloat & a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );

	l_result = a_other + *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator+=(double a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );
	VPFloat l_other(a_other);

	l_result = l_other + *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator-=(const VPFloat & a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );

	l_result = a_other - *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator-=(double a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );
	VPFloat l_other(a_other);

	l_result = l_other - *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator*=(const VPFloat & a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );

	l_result = a_other * *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator*=(double a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );
	VPFloat l_other(a_other);

	l_result = l_other * *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator/=(const VPFloat & a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );

	l_result = a_other / *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator/=(double a_other) {
	VPFloat l_result( this->m_environment.es, this->m_environment.bis, this->m_environment.stride );
	VPFloat l_other(a_other);

	l_result = l_other / *this;

	*this = l_result;

	return *this;
}

VPFloat & VPFloat::operator=(const VPFloat & a_other) {
	vassign(a_other.m_data, a_other.m_environment, this->m_data, this->m_environment);
	return *this;
}

VPFloat & VPFloat::operator=(double a_other) {
	vassign((void *)&a_other, VPFLOAT_EVP_DOUBLE, this->m_data, this->m_environment);
	return *this;
}

VPFloat VPFloat::pow2(int n, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride ) {
	VPFloat l_result(a_exponent_size, a_bis, a_stride);

	::pow2(n, l_result.m_data, l_result.m_environment);

	return l_result;
}

/*******************************************************************************************************************
 * VPFloatArray Class
 ******************************************************************************************************************/

VPFloatArray::VPFloatArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, int a_nb_elements) :
	VPFloat(a_exponent_size, a_bis, a_stride),
	m_nb_elements(a_nb_elements)
{   
	this->m_data = malloc(VPFLOAT_SIZEOF(this->m_environment) * this->m_nb_elements);

	memset(this->m_data, 0, a_nb_elements * VPFLOAT_SIZEOF(this->m_environment));

	this->m_release_m_data_on_destruction = true;
}

VPFloatArray::VPFloatArray(VPFloatArray & a_other):
	VPFloat(a_other.m_environment.es, a_other.m_environment.bis, a_other.m_environment.stride),
	m_nb_elements(a_other.m_nb_elements)
{
	this->m_data = a_other.m_data;

	this->m_release_m_data_on_destruction = false;
}

VPFloatArray::VPFloatArray(VPFloatArray * a_other):
	VPFloat(a_other->m_environment.es, a_other->m_environment.bis, a_other->m_environment.stride),
	m_nb_elements(a_other->m_nb_elements)
{
	this->m_data = a_other->m_data;

	this->m_release_m_data_on_destruction = false;
}

VPFloatArray::VPFloatArray(double * a_other, int a_nb_elements):
	VPFloat(VPFLOAT_EVP_DOUBLE.es, VPFLOAT_EVP_DOUBLE.bis, VPFLOAT_EVP_DOUBLE.stride)
{
	this->m_data = a_other;
	this->m_nb_elements = a_nb_elements;

	this->m_release_m_data_on_destruction = false;
}

VPFloatArray::VPFloatArray(float * a_other, int a_nb_elements):
	VPFloat(VPFLOAT_EVP_FLOAT.es, VPFLOAT_EVP_FLOAT.bis, VPFLOAT_EVP_FLOAT.stride)
{
	this->m_data = a_other;
	this->m_nb_elements = a_nb_elements;

	this->m_release_m_data_on_destruction = false;
}

VPFloatArray::~VPFloatArray()
{
	if ( this->m_release_m_data_on_destruction ) {
		free(this->m_data);
		this->m_data = NULL;
	}
}

void * VPFloatArray::operator+(int nb_elements) const {
	return vat(this->m_data, this->m_environment, nb_elements);
}

VPFloat VPFloatArray::operator[](int nb_elements) const {
	return VPFloat(vat(this->m_data, this->m_environment, nb_elements), this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
}

/*******************************************************************************************************************
 * VPFloatComputingEnvironment Class
 ******************************************************************************************************************/
vpfloat_ec_t VPFloatComputingEnvironment::m_environment = {0};
vpfloat_evp_t VPFloatComputingEnvironment::m_temporary_var_environment = VPFLOAT_EVP_MAX;

void VPFloatArray::printAsVector(char * a_label) const {
	vprint_vector(a_label, this->m_nb_elements, this->m_data, this->m_environment);
}
        
