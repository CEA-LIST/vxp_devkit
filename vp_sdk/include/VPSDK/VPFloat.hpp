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

#ifndef __VPFLOAT_HPP__
#define __VPFLOAT_HPP__

#include "VRPSDK/asm/vpfloat.h"
#include <iostream>

namespace VPFloatPackage {

    enum VPFloatRoundingMode { VP_RNE, VP_RTZ, VP_RDN, VP_RUP, VP_RMM };

    /*******************************************************************************************************************
     * VPFloatComputingEnvironment
     * This class is used to manage the configuration environment for vpfloat arithmetic operations.
     * There is only one Computing environment supported for all vpfloat arthmetic operations. If user need to change 
     * this environment, the modification must be done before each operation requiring this modification.
     ******************************************************************************************************************/
    class VPFloatComputingEnvironment {
        private:
            static vpfloat_ec_t m_environment;
            static vpfloat_evp_t m_temporary_var_environment;

        public:
            /*
             * Function used ,mainly by VPFloat arthmetic operators, to retrieve the arithmetic configuration environment.
             */
            static vpfloat_ec_t get_environment() { return VPFloatComputingEnvironment::m_environment; }

            /*
             * Change the precision used by the next vpfloat  arithmetic operations.
             */
            static void set_precision(uint16_t a_precision) {
                VPFloatComputingEnvironment::m_environment.wp = a_precision;
            }

            /*
             * Change the rounding mode used by the next vpfloat  arithmetic operations.
             * Possible rounding mode are:
             *   - VPFloatRoundingMode::VP_RNE
             *   - VPFloatRoundingMode::VP_RTZ
             *   - VPFloatRoundingMode::VP_RDN
             *   - VPFloatRoundingMode::VP_RUP
             *   - VPFloatRoundingMode::VP_RMM
             */
            static void set_rounding_mode(VPFloatRoundingMode a_rounding_mode);

            static uint8_t get_rounding_mode();

            /*
             * Returns the current configured precision.
             */
            static uint16_t get_precision(){
                return VPFloatComputingEnvironment::m_environment.wp;
            }

            static void set_tempory_var_environment(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride) {
                VPFloatComputingEnvironment::m_temporary_var_environment.es = a_exponent_size;
                VPFloatComputingEnvironment::m_temporary_var_environment.bis = a_bis;
                VPFloatComputingEnvironment::m_temporary_var_environment.stride = a_stride;
            }

            static vpfloat_evp_t get_temporary_var_environment() {
                return VPFloatComputingEnvironment::m_temporary_var_environment;
            }
    };

    /*******************************************************************************************************************
     * VPFloat Class
     * This class is used to manipulate vpfloat numbers. A VPFloat object is composed of two parts:
     *   - a pointer to the memory where relies the memory representation of vpfloat number
     *   - a structure defining the configuration environment to use to manipulate the vpfloat number.
     * When this code is used on a non VRP plateform, vpfloat number are mapped to double number.
     ******************************************************************************************************************/
    class VPFloat {
        public:
            /*
             * Avoid to have a default VPFloat construction to ensure the used of
             * an idedntified vpfloat environment.
             */
            VPFloat();

            VPFloat(const VPFloat & a_other);

            /*
             * This constructor allocate memory to store a vpfloat number that match  the environment defined by the parameters
             * given to the function:
             *   - exponent_size
             *   - bis
             *   - stride
             * VPFloat environment is configured accourding these values too.
             */
            VPFloat(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride );

            /*
             * This constructor do not allocate memory. It use the address given as parameters to point the vpfloat number data.
             * The VPFloat environment is configured by using the function parameters
             * given to the function:
             *   - exponent_size
             *   - bis
             *   - stride
             */            
            VPFloat( void * a_data, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, bool a_release_memory_on_destruction = false );

            /*
             * This constructor allocate memory to store a vpfloat number that match VPFLOAT_EVP_DOUBLE (defined in 
             * asm/vpfloat.h).
             * the double value given as parameter is copied into the newly allocated memory.
             * VPFloat environment is set to VPFLOAT_EVP_DOUBLE.
             */
            VPFloat( double a_value );
            
            #ifdef ENABLE_EIGEN_INTERFACE
            
            /*
             * This constructor allocate memory to store a vpfloat number that match VPFLOAT_EVP_DOUBLE (defined in 
             * asm/vpfloat.h).
             * the int value given as parameter converted into a double and then is copied into the newly allocated memory.
             * VPFloat environment is set to VPFLOAT_EVP_DOUBLE.
             */
            VPFloat( int a_value ) : VPFloat(double(a_value)) {}
           
            /*
             * This constructor allocate memory to store a vpfloat number that match VPFLOAT_EVP_DOUBLE (defined in 
             * asm/vpfloat.h).
             * the long int value given as parameter converted into a double and then is copied into the newly allocated memory.
             * VPFloat environment is set to VPFLOAT_EVP_DOUBLE.
             */
            VPFloat( long int a_value ) : VPFloat(double(a_value)) {}
            
            /*
             * This constructor allocate memory to store a vpfloat number that match VPFLOAT_EVP_DOUBLE (defined in 
             * asm/vpfloat.h).
             * the long int value given as parameter converted into a double and then is copied into the newly allocated memory.
             * VPFloat environment is set to VPFLOAT_EVP_DOUBLE.
             */
            VPFloat( long unsigned int a_value ) : VPFloat(double(a_value)) {}

						#endif // ENABLE_EIGEN_INTERFACE

            /*
             * This constructor allocate memory to store a vpfloat number that match VPFLOAT_EVP_FLOAT (defined in 
             * asm/vpfloat.h).
             * the float value given as parameter is copied into the newly allocated memory.
             * VPFloat environment is set to VPFLOAT_EVP_FLOAT.
             */
            VPFloat( float a_value );

            /*
             * This constructor allocate memory to store a vpfloat number that match the environment defined by the parameters
             * given to the function:
             *   - exponent_size
             *   - bis
             *   - stride
             * VPFloat environment is configured according these values too.
             * The value of the VPFloat is built from mantissa and exponent parameters.
             */
            VPFloat( uint64_t a_mantissa, uint64_t a_exponent, bool a_sign, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride );

            /*
             * Destructor for VPFloat objects
             */
            ~VPFloat();

            /*
             * Function used to retrieve the Vpfloat environment linked to the current VPFloat object.
             */
            vpfloat_evp_t getEnvironment() const;

            /*
             * Function used to retrieve the memory buffer where the VPFloat is stored. 
             * Needed by VLBAS.
             */
            void * getData() const;

            /*
             * This operation add the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator+=( const VPFloat & a_other);
            
            /*
             * This operation add the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator+=( double a_other);
            
            /*
             * This operation substract the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator-=( const VPFloat & a_other);
            
            /*
             * This operation substract the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator-=( const double a_other);

						/*
             * This operation multiplies the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator*=( const VPFloat & a_other);
            
					  /*
             * This operation multiplies the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator*=( const double a_other);

					  /*
             * This operation divide the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator/=( const VPFloat & a_other);
            
					  /*
             * This operation divide the value contained in the other parameter to the value of the current VPFloat and
             * store the result of the operation in the current VPFloat.
             * The environment of the result is the same than the one of the current VPFloat.
             */
            VPFloat & operator/=( const double a_other);

            /*
             * Store the value from the VPFloat number given as parameter into the current VPFloat. The value is
             * converted to the current VPFloat environment.
             */
            VPFloat & operator=( const VPFloat & a_other);

            /*
             * Store the value from the double number given as parameter into the current VPFloat. The value is
             * converted to the current VPFloat environment.
             */
            VPFloat & operator=( double a_other);

            /*
             * Function returning a double approximation of the current VPFloat number
             */
            double getdouble() const;

            /*
             * Function returning a float approximation of the current VPFloat number
             */
            float getfloat() const;

            /*
             * Function returning current VPFloat exponent value
             */
            int64_t exponent() const;

            /*
             * Function setting the current VPFloat exponent value
             */
            void exponent(const int64_t a_exponent);

            /*
             * Function returning current VPFloat mantissa chunk of index chunk_index
             */
            uint64_t mantissaChunk(const uint16_t a_chunk_index) const;

            /*
             * Function setting current VPFloat mantissa chunk of index chunk_index
             */
            void mantissaChunk(const uint16_t a_chunk_index, const uint64_t a_chunk_value);

            bool isNaN() const;
            
            bool isInf() const;

            /*******
             * Convertion functions
             ******/	
            /*
             * Maybe use explicit keyword here so that
             * 
             * VPFloat a(1.);
             * double b = a; 
             * 
             * do not work. The conversion is still possible but has to be stated explicitly:
             * double b = static_cast<double>(a);
             * 
             */
						
            /*
             * Returns the convertion of the current VPFloat number into double
             */
            operator double() const  {
                // std::cout << "conversion to double" << std::endl;   /* FIXME: add trace on debug */
                return this->getdouble();
            }

            /*
             * Returns the convertion of the current VPFloat number into float
             */
            operator float() const {
                // std::cout << "conversion to float" << std::endl;  /* FIXME: add trace on debug */
                return this->getfloat();
            }

            /*******
             * Friend operators
             ******/
            /*
             * Return the absolute value of a VPFloat scalar.
             * The result of the operator is a new VPFloat configured with same precision than the input argument
             */
            friend VPFloat abs(const VPFloat & a_x);

            /*
             * Addition operator between two VPFloat.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator+( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Addition operator between a VPFloat and a double.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator+( const VPFloat & a_lhs, double a_rhs );

            /*
             * Substraction operator between two VPFloat.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator-( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Substraction operator between a VPFloat and a double.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator-( const VPFloat & a_lhs, double a_rhs );

            /*
             * Substraction operator between a double and a VPFloat.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
	        friend VPFloat operator-(double a_lhs, const VPFloat & a_rhs);

            /*
             * Multiplication operator between two VPFloat.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator*( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Multiplication operator between a VPFloat and a double.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator*( const VPFloat & a_lhs, double a_rhs );

            /*
             * Division operator between two VPFloat.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator/( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Division operator between a VPFloat and a double.
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is a new VPFloat configured with the maximum precision (VPFLOAT_EVP_MAX)
             */
            friend VPFloat operator/( const VPFloat & a_lhs, double a_rhs );

            /*
             * Sign invertion operator.
             * The result is a new VPFloat configured with the same environment than the current VPFloat.
             */
            friend VPFloat operator-( const VPFloat & a_rhs);

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is greater than a_rhs
             */
            friend bool operator>( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is greater than or equal to a_rhs
             */
            friend bool operator>=( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Comparison of VPFloat scalar to a double
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is smaller than a_rhs
             */
            friend bool operator<( const VPFloat & a_lhs, double a_rhs );

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is smaller than a_rhs
             */
            friend bool operator<( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is smaller than or equal to a_rhs
             */
            friend bool operator<=( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is equal to a_rhs
             */
            friend bool operator==( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Comparison of two VPFloat scalar
             * This operator use to current configuration set into VPFloatComputingEnvironment.
             * The result of the operator is true if a_lhs is not equal to a_rhs
             */
            friend bool operator!=( const VPFloat & a_lhs, const VPFloat & a_rhs );

            /*
             * Function used to print VPFloat number
             */
            friend std::ostream& operator << (std::ostream & stream, const VPFloat & a_vpfloat);

            /*******
             * Friend components
             *******/
            friend class VPFloatOperation;

            /*******
             * Static functions
             *******/
            static VPFloat pow2(int n, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride );

        protected:
            /*
             * Configuration environment for the current VPFloat
             */
            vpfloat_evp_t m_environment;

            /*
             * Pointer to the memory representation of the current VPFloat.
             */
            void * m_data;

            /*
             * Flag set when current object is loaded with a double number.
             */
            bool m_double_data;

            /*
             * Flag set when current object is loaded with a float number.
             */
            bool m_float_data;

            /*
             * Flag set when m_data was allocated by a class function
             */
            bool m_release_m_data_on_destruction;
    };

    /*******************************************************************************************************************
     * VPFloatArray Class
     * This class is used to manage array of vpfloats.
     ******************************************************************************************************************/
    class VPFloatArray: public VPFloat {
        public:
            /*
             * Constructor allocating memory of vpfloat number matching the configuration provided throuh parameters:
             *   - exponent_size
             *   - bis
             *   - stride
             * The array will contain nb_elements contiguous vpfloat.
             */
            VPFloatArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, int a_nb_elements);

            /*
             * Constructor used to copy one array to a new one.
             * Data are copied.
             */
            VPFloatArray(VPFloatArray & a_other);

            /*
             * Constructor used to copy one array to a new one.
             * Data are copied.
             */
            VPFloatArray(VPFloatArray * a_other);

            /*
             * Constructor used to build a VPFloatArray upon an existing double array.
             */
            VPFloatArray(double * a_data, int a_nb_elements);

            /*
             * Constructor used to build a VPFloatArray upon an existing float array.
             */
            VPFloatArray(float * a_data, int a_nb_elements);

            /*
             * Release memory allocated for the vpfloat array
             */
            ~VPFloatArray();

            /*
             * This operator returns the address of the element at index of the array.
             */
            void * operator+(int nb_elements) const;

            /*
             * This operator return a VPFloat object that point to the data of the element at index of the array.
             * If modification are made to the returned VPFloat, the data in the array are modified.
             */
            VPFloat operator [](int index) const;

            /*
             * Function returning the number of element stored in the array
             */
            int nbElements() const{
                return this->m_nb_elements;
            }

            /*
             * Function displaying VPFloatArray content as a vector
             */
            void printAsVector(char * a_label) const;

        private:
            /*
             * Number of elements in the array.
             */
            int m_nb_elements;
    };


    /* Function returning abolute value of a VPFloat */
    VPFloat abs(const VPFloat & a_x);

}

#endif /* __VPFLOAT_HPP__ */
