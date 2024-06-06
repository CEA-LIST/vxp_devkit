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
#include <iomanip>
#include <cmath>
#include <cstddef>
#include <math.h>
#include <mpfr.h>

using namespace VPFloatPackage;

/**
 * < when writing a vp number into a char array, following the notation "0b1.[significand]p[exponent in decimal]", 
 * max number of chars is 
 * 1 (sign) + 4 (0b1.) + P_CHUNK_LEN*P_CHUNK_MAX_COUNT (significand) + 2 (exponent symbol + its sign) + 18 (58-bits signed exponent in decimal) + 1 (end of string)
 **/
const uint16_t MAX_MPFR_STRING_LEN = 26+P_CHUNK_LEN*P_CHUNK_MAX_COUNT;
const uint16_t MPFR_STRING_EXPONENT_LEN = 20;


/*******************************************************************************************************************
 * VPFloatComputingEnvironment
 ******************************************************************************************************************/
void VPFloatComputingEnvironment::set_rounding_mode(VPFloatRoundingMode a_rounding_mode) {
	switch(a_rounding_mode) {
		case VP_RNE:
			mpfr_set_default_rounding_mode(MPFR_RNDN);
			break;		
		case VP_RTZ:
			mpfr_set_default_rounding_mode(MPFR_RNDZ);
			break;		
		case VP_RDN:
			mpfr_set_default_rounding_mode(MPFR_RNDD);
			break;		
		case VP_RUP:
			mpfr_set_default_rounding_mode(MPFR_RNDU);
			break;		
		case VP_RMM:
			mpfr_set_default_rounding_mode(MPFR_RNDN);
			std::cout << "VPFloat roundingmode VP_RMM converted to MPFR_RNDN" << std::endl << std::flush;
			break;
	}
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

	this->m_data = (mpfr_t *)malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;

	mpfr_init2(*((mpfr_t *)(this->m_data)), this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpfr_set_d(*((mpfr_t *)(this->m_data)), 0.0, mpfr_get_default_rounding_mode());
	mpfr_set_exp(*((mpfr_t *)(this->m_data)), (mpfr_exp_t)0);
}

VPFloat::VPFloat(const VPFloat & a_other) :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_other.m_environment.es;
	this->m_environment.bis = a_other.m_environment.bis;
	this->m_environment.stride = a_other.m_environment.stride;

	this->m_data = (mpfr_t *)malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;

    mpfr_init2(*((mpfr_t *)(this->m_data)),mpfr_get_prec(*((mpfr_t *)(a_other.m_data))));
    mpfr_set  (*((mpfr_t *)(this->m_data)),*((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());
}

VPFloat::VPFloat(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride) :
	m_double_data(false),
	m_float_data(false)
{
	this->m_environment.es = a_exponent_size;
	this->m_environment.bis = a_bis;
	this->m_environment.stride = a_stride;

	this->m_data = (mpfr_t *)malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;

	mpfr_init2(*((mpfr_t *)(this->m_data)), this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpfr_set_d(*((mpfr_t *)(this->m_data)), 0.0, mpfr_get_default_rounding_mode());
	mpfr_set_exp(*((mpfr_t *)(this->m_data)), (mpfr_exp_t)0);
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

	this->m_data = malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;

	mpfr_init2(*((mpfr_t *)(this->m_data)), this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpfr_set_d(*((mpfr_t *)(this->m_data)), 1.0, mpfr_get_default_rounding_mode());
	if ( mpfr_set_exp(*((mpfr_t *)(this->m_data)), a_exponent + 1) != 0 ) {
		std::cout << "Fail setting exponent value in constructor." << std::endl;
	}	

	this->mantissaChunk(0, a_mantissa);

}

VPFloat::VPFloat( double a_value ) :
	m_double_data(true),
	m_float_data(false)
{
	this->m_environment.es = VPFLOAT_EVP_DOUBLE.es;
	this->m_environment.bis = VPFLOAT_EVP_DOUBLE.bis;
	this->m_environment.stride = VPFLOAT_EVP_DOUBLE.stride;

	this->m_data = malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;

	mpfr_init2(*((mpfr_t *)(this->m_data)), this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpfr_set_d(*((mpfr_t *)(this->m_data)), a_value, mpfr_get_default_rounding_mode());

}

VPFloat::VPFloat( float a_value ) :
	m_double_data(false),
	m_float_data(true)
{
	this->m_environment.es = VPFLOAT_EVP_FLOAT.es;
	this->m_environment.bis = VPFLOAT_EVP_FLOAT.bis;
	this->m_environment.stride = VPFLOAT_EVP_FLOAT.stride;

	this->m_data = malloc(sizeof(mpfr_t));
	this->m_release_m_data_on_destruction = true;
	mpfr_init2(*((mpfr_t *)(this->m_data)), this->m_environment.bis - this->m_environment.es - 1 + 1);
	mpfr_set_d(*((mpfr_t *)(this->m_data)), (double)a_value, mpfr_get_default_rounding_mode());
}

VPFloat::~VPFloat() {
	if ( this->m_release_m_data_on_destruction && this->m_data != NULL ) {
		mpfr_clear(*(mpfr_t *)(this->m_data));
		free(this->m_data);
	}
}

vpfloat_evp_t VPFloat::getEnvironment() const {
	return this->m_environment;
}

void * VPFloat::getData() const {
	return this->m_data;
}

int64_t VPFloat::exponent() const
{
	return (int64_t) (mpfr_get_exp(*((mpfr_t *)(this->m_data))) - 1);
}

double VPFloat::getdouble() const
{
	return mpfr_get_d(*((mpfr_t *)(this->m_data)), mpfr_get_default_rounding_mode());
}

float VPFloat::getfloat() const
{
	return (float)mpfr_get_d(*((mpfr_t *)(this->m_data)), mpfr_get_default_rounding_mode());
}

void VPFloat::exponent(const int64_t a_exponent)
{
	int l_rc = mpfr_set_exp(*((mpfr_t *)(this->m_data)), a_exponent + 1);

	if ( l_rc != 0 ) {
		std::cout << "Fail setting VPFloat exponent to " << a_exponent << " rc is " << l_rc << std::endl;
	}
}

uint64_t VPFloat::mantissaChunk(const uint16_t a_chunk_index) const {
	mpfr_exp_t l_exponent = 0;
	uint64_t l_chunk_value = 0x0000000000000000;
	int l_bit_offset_for_chunk_end, l_bit_offset_for_chunk, l_max_offset;

	/*
	 * Use MPFR function to retrieve the string containing the 
	 *  binary representation of the number
	 */
	char * l_str_buffer = mpfr_get_str(NULL, &l_exponent, 2, 0, (*(mpfr_t *)(this->m_data)), mpfr_get_default_rounding_mode());

	if ( l_str_buffer == NULL ) {
		std::cout << "Fail retrieving string representing VPFloat number." << std::endl;
		return 0;
	}

	/*
	 * Compute the offeset of the data representing the requested 
	 * chunk (skip first bit to)
	 */
	int l_string_offset = 1; // Skip first bit

	if ( l_str_buffer[0] == '-' ) {
		l_string_offset++; // Skip sign bit
	}

	l_bit_offset_for_chunk_end = l_string_offset + ( a_chunk_index + 1) * P_CHUNK_LEN;
	l_bit_offset_for_chunk = l_string_offset + a_chunk_index * P_CHUNK_LEN;
	l_max_offset = strlen(l_str_buffer);

	/*
	 * Loop in chunk string representation to build chunk numeric value
	 */
	for ( int l_current_offset = l_bit_offset_for_chunk; ( l_current_offset < l_bit_offset_for_chunk_end ) && ( l_current_offset <= l_max_offset ) ; l_current_offset++ ) {
		int l_current_bit_index = l_bit_offset_for_chunk_end - l_current_offset - 1;
		if ( l_str_buffer[ l_current_offset ] == '1' ) {
			l_chunk_value = l_chunk_value | ( 1ULL << l_current_bit_index );
		} else if ( l_str_buffer[ l_current_offset ] == '\0' ) { 
			break;
		}

		l_bit_offset_for_chunk++;
	}

	
	/*
	 * Release initial string buffer used by mpfr_get_str.
	 */
	mpfr_free_str(l_str_buffer);

	/*
	 * Retrieve the uint64_t object with the nartive bitset function.
	 */
	return l_chunk_value;
}

void VPFloat::mantissaChunk(const uint16_t a_chunk_index, const uint64_t a_chunk_value) {
	
	mpfr_exp_t l_exponent = 0;
	char * l_current_str_value = (char *)malloc(sizeof(char) * MAX_MPFR_STRING_LEN);
	char * l_new_str_value     = (char *)malloc(sizeof(char) * MAX_MPFR_STRING_LEN);
	int l_string_offset = 4; // shift offset for "0b1." that will be prepend to the string
	int l_bit_offset_for_chunk_end, l_bit_offset_for_chunk;
	int l_initial_str_length = 0;

	/*
	 * Use MPFR function to retrieve the string containing the 
	 * binary representation of the number
	 */
	if ( mpfr_get_str(l_current_str_value, &l_exponent, 2, 0, (*(mpfr_t *)(this->m_data)), mpfr_get_default_rounding_mode()) == NULL ) {
		std::cout << "Fail retrieving string representing VPFloat number. Can't set matissa chunk " << a_chunk_index << std::endl;
		return;	
	}

 	l_initial_str_length = strlen(l_current_str_value);

	/*
	 * Shift string offset since we need to consider the sign
	 */
	if( l_current_str_value[0] == '-' ) {
		l_string_offset += 1;
	}

	/*
	 * Compute the offeset of the data representing the requested 
	 * chunk (skip first bit to)
	 */
	l_bit_offset_for_chunk_end = l_string_offset + ( a_chunk_index + 1) * P_CHUNK_LEN;
	l_bit_offset_for_chunk = l_string_offset + a_chunk_index * P_CHUNK_LEN;

	/*
	 * Allocate a buffer large enough to receive additional chunk infirmations
	 */
	if( l_current_str_value[0] == '-' ) {
		sprintf(l_new_str_value, "-0b1.%s", l_current_str_value + 1);
	} else {
		sprintf(l_new_str_value,  "0b1.%s", l_current_str_value + 1);
	}

	/*
	 * Release initial string buffer used by mpfr_get_str since data were copied to the
	 * new buffer.
	 */
	mpfr_free_str(l_current_str_value);

	/*
	 * Replace desired chunk in intial string representation of the number.
	 */
	for ( int l_current_offset = l_bit_offset_for_chunk; l_current_offset < l_bit_offset_for_chunk_end; l_current_offset++ ) {
		int l_current_bit_index = ( l_bit_offset_for_chunk_end - l_current_offset ) - 1;
		uint64_t l_bit_value = ( a_chunk_value >> l_current_bit_index ) & 1ULL ;
		if ( l_bit_value == 1 ) {
			l_new_str_value[ l_current_offset ] = '1';
		} else {
			l_new_str_value[ l_current_offset ] = '0';
		}
	}

	/*
	 * Append the exponent value at the end of the string
	 */
	sprintf(l_new_str_value + std::max(l_initial_str_length, l_bit_offset_for_chunk_end), "p%ld", l_exponent-1);

	/*
	 * Convert the string containing the binary representation of
	 * the number into the MPFR number
	 */
	int l_rc = mpfr_set_str((*(mpfr_t *)(this->m_data)), l_new_str_value, 2, mpfr_get_default_rounding_mode());

	if ( l_rc != 0 ) {
		std::cout << "Fail setting VPFloat from its string representation." << std::endl;
	}

	free(l_new_str_value);
}

bool VPFloat::isNaN() const {
	return ( mpfr_nan_p((*(mpfr_t *)(this->m_data))) != 0 );
}

bool VPFloat::isInf() const {
	return ( mpfr_inf_p((*(mpfr_t *)(this->m_data))) != 0 );
}

namespace VPFloatPackage {

	VPFloat abs(const VPFloat & a_x) {
		VPFloat l_abs(a_x.m_environment.es, a_x.m_environment.bis, a_x.m_environment.stride);

		mpfr_abs(*((mpfr_t *)(l_abs.m_data)), *((mpfr_t *)(a_x.m_data)), mpfr_get_default_rounding_mode());

		return l_abs;
	}

	VPFloat operator+(const VPFloat & a_lhs,const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);
		
		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_add(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator+( const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_add_d(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), a_rhs, mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator-(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);
		
		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_sub(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());
		return l_result;
	}

	VPFloat operator-(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_sub_d(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), a_rhs, mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator-(double a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_t l_lhs;
		mpfr_init_set_d(l_lhs, a_lhs, mpfr_get_default_rounding_mode());
		mpfr_sub(*((mpfr_t *)(l_result.m_data)), l_lhs, *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator*(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);
		
		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_mul(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator*(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_mul_d(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), a_rhs, mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator/(const VPFloat & a_lhs, const VPFloat & a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_div(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator/(const VPFloat & a_lhs, double a_rhs) {
		vpfloat_evp_t  l_tmp_var_environment = VPFloatComputingEnvironment::get_temporary_var_environment();
		VPFloat l_result(l_tmp_var_environment.es, l_tmp_var_environment.bis, l_tmp_var_environment.stride);

		mpfr_set_default_prec(VPFloatComputingEnvironment::get_precision());

		mpfr_div_d(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_lhs.m_data)), a_rhs, mpfr_get_default_rounding_mode());

		return l_result;
	}

	VPFloat operator-(const VPFloat & a_rhs) {
		VPFloat l_result(a_rhs.m_environment.es, a_rhs.m_environment.bis, a_rhs.m_environment.stride);
		mpfr_neg(*((mpfr_t *)(l_result.m_data)), *((mpfr_t *)(a_rhs.m_data)), mpfr_get_default_rounding_mode());
		return l_result;
	}

	bool operator>( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_greater_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) > 0;
	}

	bool operator>=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_greaterequal_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) > 0;
	}

	bool operator<( const VPFloat & a_lhs, double a_rhs ) {
		return mpfr_cmp_d(*((mpfr_t *)(a_lhs.m_data)), a_rhs) < 0;
	}

	bool operator<( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_less_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) > 0;
	}

	bool operator<=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_lessequal_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) > 0;
	}

	bool operator==( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_equal_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) > 0;
	}

	bool operator!=( const VPFloat & a_lhs, const VPFloat & a_rhs ) {
		return mpfr_equal_p(*((mpfr_t *)(a_lhs.m_data)), *((mpfr_t *)(a_rhs.m_data))) == 0;
	}

	std::ostream& operator << (std::ostream & stream, const VPFloat & a_vpfloat) {
		int l_string_width = ((a_vpfloat.m_environment.bis - a_vpfloat.m_environment.es - a_vpfloat.m_environment.stride) / 4 ) + 10;
		char * l_buffer = (char *)malloc(sizeof(char) * l_string_width );

		int l_rc = mpfr_snprintf(l_buffer, l_string_width, "%*RA", (a_vpfloat.m_environment.bis - a_vpfloat.m_environment.es - a_vpfloat.m_environment.stride)/4 , *((mpfr_t *)(a_vpfloat.m_data)));
		
		if ( l_rc < 0 ) {
			stream << "Fail calling mpfr_snprintf to display VPFloat." << std::endl;
		} else {
			stream << l_buffer << std::flush;
		}

		free(l_buffer);
			
		return stream;
	}
}

VPFloat & VPFloat::operator+=(const VPFloat & a_other) {
	mpfr_add(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), *((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator+=(double a_other) {
	mpfr_add_d(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), a_other, mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator-=(const VPFloat & a_other) {
	mpfr_sub(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), *((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator-=(double a_other) {
	mpfr_sub_d(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), a_other, mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator*=(const VPFloat & a_other) {
	mpfr_mul(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), *((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator*=(double a_other) {
	mpfr_mul_d(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), a_other, mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator/=(const VPFloat & a_other) {
	mpfr_div(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), *((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator/=(double a_other) {
	mpfr_div_d(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(this->m_data)), a_other, mpfr_get_default_rounding_mode());

	return *this;
}

VPFloat & VPFloat::operator=(const VPFloat & a_other) {
	mpfr_set(*((mpfr_t *)(this->m_data)), *((mpfr_t *)(a_other.m_data)), mpfr_get_default_rounding_mode());
	return *this;
}

VPFloat & VPFloat::operator=(double a_other) {
	mpfr_set_d(*((mpfr_t *)(this->m_data)), a_other, mpfr_get_default_rounding_mode());
	return *this;
}

VPFloat VPFloat::pow2(int n, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride ) {
	VPFloat l_result(a_exponent_size, a_bis, a_stride);

	mpfr_ui_pow_ui(*((mpfr_t *)(l_result.m_data)), 2, n, mpfr_get_default_rounding_mode());

	return l_result;
}

/*******************************************************************************************************************
 * VPFloatArray Class
 ******************************************************************************************************************/
VPFloatArray::VPFloatArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, int a_nb_elements) :
	VPFloat(a_exponent_size, a_bis, a_stride),
	m_nb_elements(a_nb_elements)
{   
	int l_index = 0;
	// Release m_data buffer allocated by VPFloat since we overwrite it to store our array
	mpfr_clear(*((mpfr_t *)this->m_data));
	free(this->m_data);
	this->m_data = malloc(sizeof(mpfr_t) * this->m_nb_elements);
	for (l_index = 0 ; l_index < a_nb_elements; l_index++) {
		mpfr_init2(((mpfr_t *)this->m_data)[l_index], a_bis - a_exponent_size - 1 + 1);
		mpfr_set_d(((mpfr_t *)this->m_data)[l_index], 0.0, mpfr_get_default_rounding_mode());
		mpfr_set_exp(((mpfr_t *)this->m_data)[l_index], (mpfr_exp_t)0);
	}

	this->m_release_m_data_on_destruction = true;
}

VPFloatArray::VPFloatArray(VPFloatArray & a_other):
	VPFloat(a_other.m_environment.es, a_other.m_environment.bis, a_other.m_environment.stride),
	m_nb_elements(a_other.m_nb_elements)
{
	// Release m_data buffer allocated by VPFloat since we overwrite it to store our array
	mpfr_clear(*((mpfr_t *)this->m_data));
	free(this->m_data);
	this->m_data = a_other.m_data;

	this->m_release_m_data_on_destruction = false;
}

VPFloatArray::VPFloatArray(VPFloatArray * a_other):
	VPFloat(a_other->m_environment.es, a_other->m_environment.bis, a_other->m_environment.stride),
	m_nb_elements(a_other->m_nb_elements)
{
	// Release m_data buffer allocated by VPFloat since we overwrite it to store our array
	mpfr_clear(*((mpfr_t *)this->m_data));
	free(this->m_data);
	this->m_data = a_other->m_data;

	this->m_release_m_data_on_destruction = false;
}


VPFloatArray::VPFloatArray(double * a_other, int a_nb_elements):
	VPFloat(VPFLOAT_EVP_DOUBLE.es, VPFLOAT_EVP_DOUBLE.bis, VPFLOAT_EVP_DOUBLE.stride),
	m_nb_elements(a_nb_elements)
{
	int l_index = 0;

	// Release m_data buffer allocated by VPFloat since we overwrite it to store our array
	mpfr_clear(*((mpfr_t *)this->m_data));
	free(this->m_data);
	this->m_data = malloc(sizeof(mpfr_t) * this->m_nb_elements);
	for (l_index = 0 ; l_index < a_nb_elements; l_index++) {
		mpfr_init2(((mpfr_t *)this->m_data)[l_index], VPFLOAT_EVP_DOUBLE.bis - VPFLOAT_EVP_DOUBLE.es - 1 + 1);
		mpfr_set_d(((mpfr_t *)this->m_data)[l_index], a_other[l_index], mpfr_get_default_rounding_mode());
	}

	this->m_release_m_data_on_destruction = true;
}


VPFloatArray::VPFloatArray(float * a_other, int a_nb_elements):
	VPFloat(VPFLOAT_EVP_FLOAT.es, VPFLOAT_EVP_FLOAT.bis, VPFLOAT_EVP_FLOAT.stride),
	m_nb_elements(a_nb_elements)
{
	int l_index = 0;

	// Release m_data buffer allocated by VPFloat since we overwrite it to store our array
	mpfr_clear(*((mpfr_t *)this->m_data));
	free(this->m_data);
	this->m_data = malloc(sizeof(mpfr_t) * this->m_nb_elements);
	for (l_index = 0 ; l_index < a_nb_elements; l_index++) {
		mpfr_init2(((mpfr_t *)this->m_data)[l_index], VPFLOAT_EVP_FLOAT.bis - VPFLOAT_EVP_FLOAT.es - 1 + 1);
		mpfr_set_d(((mpfr_t *)this->m_data)[l_index], a_other[l_index], mpfr_get_default_rounding_mode());
	}

	this->m_release_m_data_on_destruction = true;
}

VPFloatArray::~VPFloatArray()
{
	if ( this->m_release_m_data_on_destruction ) {
		int l_index = 0;
		for (l_index = 0 ; l_index < this->m_nb_elements; l_index++) {
			mpfr_clear(((mpfr_t *)this->m_data)[l_index]);
		}
		free(this->m_data);
		this->m_data = NULL;
	}
}

void * VPFloatArray::operator+(int nb_elements) const {
	mpfr_t * l_mpfr_address = (mpfr_t *)((unsigned long long)(this->m_data) + (nb_elements * sizeof(mpfr_t)));
	return (void *)l_mpfr_address;
}

VPFloat VPFloatArray::operator[](int nb_elements) const {
	mpfr_t * l_mpfr_address = (mpfr_t *)((unsigned long long)(this->m_data) + (nb_elements * sizeof(mpfr_t)));
	return VPFloat((void *)l_mpfr_address, this->m_environment.es, this->m_environment.bis, this->m_environment.stride);
}

void VPFloatArray::printAsVector(char * a_label) const {
	for (int l_index = 0 ; l_index < this->m_nb_elements; l_index++) {
		mpfr_t * l_mpfr_address = (mpfr_t *)((unsigned long long)(this->m_data) + (l_index * sizeof(mpfr_t)));
		double l_value = mpfr_get_d(*l_mpfr_address, mpfr_get_default_rounding_mode());
		printf("%s[%d] = %f\n", a_label, l_index, l_value);		
	}
}

/*******************************************************************************************************************
 * VPFloatComputingEnvironment Class
 ******************************************************************************************************************/
vpfloat_ec_t VPFloatComputingEnvironment::m_environment = {0};
vpfloat_evp_t VPFloatComputingEnvironment::m_temporary_var_environment = VPFLOAT_EVP_MAX;
