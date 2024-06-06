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
 * Authors       : Alexandre Hoffmann
 * Creation Date : August, 2023
 * Description   : 
 **/

#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <VPSDK/VPFloat.hpp>

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */

using Eigen::MatrixXd;
using Eigen::Matrix;

typedef Matrix<VPFloatPackage::VPFloat, 10, 10> MyMatrix;


int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;

  MyMatrix mm;
  for (int i = 0 ; i < 10; i++) {
    mm(i,i) = VPFloatPackage::VPFloat(1.0);
  }

  std::cout << mm << std::endl;

  return 0;
}