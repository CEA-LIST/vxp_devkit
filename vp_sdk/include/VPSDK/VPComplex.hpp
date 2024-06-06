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

#ifndef __VPCOMPLEX_HPP__
#define __VPCOMPLEX_HPP__

#include "VPSDK/VPFloat.hpp"
#include <iostream>

namespace VPFloatPackage {

    class VPComplex : public VPFloat {

        public:
            VPComplex();

            /*
            * Constructor allocating memory for a vpcomplex matching the configuration provided through parameters:
            *   - exponent_size
            *   - bis
            *   - stride
            * The complex will contain 2 contiguous vpfloat.
            */
            VPComplex(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride);

            /*
            * Constructor used to copy one complex to a new one.
            * Data are copied.
            */
            VPComplex(const VPComplex & a_other);

            /*
             *
             */
            VPComplex(void * a_data, vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, bool a_release_memory_on_destruction = false );

            /*
            * Release memory allocated for the vpfloat array
            */
            ~VPComplex();

            /*
             * Function used to retrieve the VPComplex environment linked to the current VPComplex object.
             */
            vpfloat_evp_t getEnvironment() const { return m_environment; }

            /*
             * Function used to retrieve the memory buffer where the VPComplex is stored. 
             * Needed by VLBAS.
             */
            void * getData() const { return m_data; }

            /*
             * Function returning the real component of Complex number
             */
            const VPFloat real() const;

            /*
             * Function returning the imaginary component of Complex number
             */
            const VPFloat imag() const;

            /*
             * Function used to set the read component of Complex number
             */
            void real(double a_x);

            /*
             * Function used to set the imaginary component of Complex number
             */
            void imag(double a_y);

            /*
             * Function returning the conjugate value of current Complex number
             */
            VPComplex conjugate() const;

            /*
             * Function returning the norm value of current Complex number
             */
            VPFloat norm() const;

            VPComplex& operator =  (const int x);
            VPComplex& operator =  (const VPComplex &a_other);
            VPComplex& operator += (const VPComplex &a_other) { this->add(a_other); return *this; }
            VPComplex& operator -= (const VPComplex &a_other) { this->sub(a_other); return *this; }
            VPComplex& operator *= (const VPComplex &a_other) { this->mul(a_other); return *this; }
            VPComplex& operator /= (const VPComplex &a_other) { this->div(a_other); return *this; }

            friend bool operator==( const VPComplex & a_lhs, const VPComplex & a_rhs );
            friend bool operator!=( const VPComplex & a_lhs, const VPComplex & a_rhs );
            friend VPComplex operator+ (const VPComplex& a_lhs, const VPComplex& a_rhs);
            friend VPComplex operator- (const VPComplex& a_lhs, const VPComplex& a_rhs);
            friend VPComplex operator* (const VPComplex& a_lhs, const VPComplex& a_rhs);
            friend VPComplex operator* (const VPComplex& a_lhs, const VPFloat& a_rhs);
            friend VPComplex operator/ (const VPComplex& a_lhs, const VPComplex& a_rhs);
            friend std::ostream& operator << (std::ostream & stream, const VPComplex & a_complex);

        private:
            /*
             *
             */
            void set_real(void * a_vpfloat_data);

            /*
             *
             */
            void set_imag(void * a_vpfloat_data);

            /*
             *  
             */
            void add(const VPComplex & a_other);

            /*
             *  
             */
            void sub(const VPComplex & a_other);

            /*
             *  
             */
            void mul(const VPComplex & a_other);

            /*
             *  
             */
            void div(const VPComplex & a_other);

            /*
             *
             */
            void init_from_environment();

            /*
             * Configuration environment for the current VPComplex
             */
            vpfloat_evp_t m_environment;

            /*
             * Pointer to the memory representation of the current VPComplex.
             */
            void * m_data;

            /*
             * Flag set when m_data was allocated by a class function
             */
            bool m_release_m_data_on_destruction;            
    };


    /*******************************************************************************************************************
     * VPFloatArray Class
     * This class is used to manage array of vpfloats.
     ******************************************************************************************************************/
    class VPComplexArray {
        public:
            VPComplexArray(vpfloat_es_t a_exponent_size, vpfloat_prec_t a_bis, vpfloat_off_t a_stride, int a_nb_elements);

            VPComplexArray(VPComplexArray & a_other);

            VPComplexArray(VPComplexArray * a_other);
            
            VPComplexArray(float * a_data, int a_nb_elements);

            VPComplexArray(double * a_data, int a_nb_elements);
            
            ~VPComplexArray();

            /*
             * Function used to retrieve the VPComplexArray environment linked to the current VPComplexArray object.
             */
            vpfloat_evp_t getEnvironment() const { return m_environment; }

            /*
             * Function used to retrieve the memory buffer where the VPComplexArray is stored. 
             * Needed by VLBAS.
             */
            void * getData() const { return m_data; }

            /*
             * This operator returns the address of the element at index of the array.
             */
            void * operator+(int nb_elements) const;

            VPComplex operator [](int index) const;

            int nbElements() const{ return this->m_nb_elements; }
            
        private:
            /*
             * Configuration environment for the current VPComplex
             */
            vpfloat_evp_t m_environment;

            /*
             * Pointer to the memory representation of the current VPComplex.
             */
            void * m_data;

            /*
             * Flag set when m_data was allocated by a class function
             */
            bool m_release_m_data_on_destruction;   

            /*
             * Number of elements in the array.
             */
            int m_nb_elements;

            uint64_t scalarSize() const;

    };

};


#endif /* __VPCOMPLEX_HPP__ */