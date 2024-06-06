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

#include <cstdlib>
#include <iostream>
#include <complex>
#include <VPSDK/VPFloat.hpp>

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */

namespace VPFloatPackage {

    class VPComplex {
        public:
            VPComplex(const VPFloat & real, const VPFloat & imag) : m_real(real), m_imag(imag) {}
            ~VPComplex() {}

            const VPFloat& real() const {return this->m_real;}
            const VPFloat& imag() const {return this->m_imag;}

            VPFloat& real() {return this->m_real;}
            VPFloat& imag() {return this->m_imag;}


        private: 
            VPFloat m_real;
            VPFloat m_imag;
    };

    const VPFloat real(const VPComplex & z) {return  z.real();}
    const VPFloat imag(const VPComplex & z) {return  z.imag();}

}


std::ostream & operator<<(std::ostream & out, const VPFloatPackage::VPComplex& z) {
    out << "(" << double(z.real()) << ", " << double(z.imag()) << ")";
    return out;
}

VPFloatPackage::VPComplex operator*(const VPFloatPackage::VPComplex & lhs, const VPFloatPackage::VPComplex & rhs) {
    return VPFloatPackage::VPComplex(lhs.real() * rhs.real() - lhs.imag() * rhs.imag(), lhs.real() * rhs.imag() - lhs.imag() * rhs.real());
}


int main(int argc, char *argv[])
{
    VPFloatPackage::VPFloatComputingEnvironment::set_precision(512);
    VPFloatPackage::VPFloatComputingEnvironment::set_tempory_var_environment(7, 53+7+1, 1);


    VPFloatPackage::VPComplex a(VPFloatPackage::VPFloat(-1.), VPFloatPackage::VPFloat(2.));
    VPFloatPackage::VPComplex b(VPFloatPackage::VPFloat(1.), VPFloatPackage::VPFloat(-2.));

    if ( double((a*b).real()) == 3 ) {
        std::cout << "SUCCESS" << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cout << "FAIL" << " (a*b).real()=" << double((a*b).real()) << std::endl;
        return EXIT_FAILURE;
    }
    
}

