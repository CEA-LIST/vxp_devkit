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

#include "VPSDK/VMath.hpp"

using namespace VPFloatPackage;

unsigned int vsqrt_segments(int p) {
  unsigned int iters=0;
  if (p<6)       iters=6;
  else if (p<20)  iters=7;
  else if (p<40)  iters=8;
  else if (p<100) iters=9;
  else if (p<180) iters=10;
  else if (p<360) iters=11;
  else if (p<400) iters=13;
  else iters=14;
  return(iters);
}

VPFloat VMath::vsqrt(const VPFloat & a) {
#ifdef DEBUG1
  std::cout<<".. entering vsqrt_recip( "<<VPFloatComputingEnvironment::get_precision()<<" , " <<a<<" ) \n";
#endif

  short myBis=VPFloatComputingEnvironment::get_precision()+a.getEnvironment().es+1;
  VPFloat xprime (a.getEnvironment().es, myBis, a.getEnvironment().stride );
  VPFloat y (a.getEnvironment().es, myBis, a.getEnvironment().stride );

  short is_odd; int c;
  int expo_a=a.exponent();

  if ((expo_a & 1)==0) { // even
    is_odd=0;
  } else {
    is_odd=1;
  }

  c=(expo_a + is_odd)/2; // ceil(exp_a/2)


  xprime=a; // copie
  xprime.exponent(expo_a-2*c);

  y=1.0;

  for (unsigned int m=0; m<vsqrt_segments(VPFloatComputingEnvironment::get_precision()); m++) {

    y=y*(VPFloat(3.0)-(xprime*y*y));

    // div y by 2
    y.exponent(y.exponent()-1);

#ifdef DEBUG1
    std::cout<<".... loop inside "<<m<<" -> " <<y<<" ) \n";
#endif
  }

  y.exponent(y.exponent()+c);

#ifdef DEBUG
    std::cout<<".... after y.exponent="<<y.exponent()<<" ) \n";
    std::cout<<".... after y="<<y<<" ) \n";
    std::cout<<".... after xprime="<<xprime<<" ) \n";
#endif

  y=xprime*y;

#ifdef DEBUG1
    std::cout<<".... after y="<<y<<" ) \n";
#endif

  return y;

}

VPFloat VMath::getMachineEpsilon()
{
	const VPFloat two(2.);
	
	VPFloat epsilon(1.);
	for (uint16_t i=0;i<VPFloatPackage::VPFloatComputingEnvironment::get_precision()-1;++i)
	{
		epsilon /= two;
	}
	return epsilon;
}

VPFloat VMath::getDummyPrecision()
{
	const VPFloat two(2.);
	
	const uint16_t weak_prec = 0.9*(VPFloatPackage::VPFloatComputingEnvironment::get_precision()-1);
	
	VPFloat epsilon(1.);
	for (uint16_t i=0;i<weak_prec-1;++i)
	{
		epsilon /= two;
	}
	return epsilon;
}
	
bool VMath::isAlmostEqual(const VPFloat& a, const VPFloat& b, const VPFloat& eps)
{
	return abs(a-b) < eps*std::min(abs(a), abs(b));
}
