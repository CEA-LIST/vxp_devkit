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

#ifndef EIGEN_VPFLOAT_SUPPORT_HPP
#define EIGEN_VPFLOAT_SUPPORT_HPP

#define ENABLE_EIGEN_INTERFACE

#include <Eigen/Core>
#include <VPSDK/VPFloat.hpp>
#include <VPSDK/VMath.hpp>

namespace Eigen 
{

template<> struct NumTraits<VPFloatPackage::VPFloat> : GenericNumTraits<VPFloatPackage::VPFloat>
{
	enum 
	{
		IsInteger = 0,
		IsSigned = 1,
		IsComplex = 0,
		RequireInitialization = 1,
		ReadCost = HugeCost,
		AddCost  = HugeCost,
		MulCost  = HugeCost
	};

	typedef VPFloatPackage::VPFloat Real;
	typedef VPFloatPackage::VPFloat NonInteger;

	static inline Real epsilon() { return VPFloatPackage::VMath::getMachineEpsilon(); }
	static inline Real dummy_precision() { return VPFloatPackage::VMath::getDummyPrecision(); }
};

}

namespace VPFloatPackage
{
	inline const VPFloatPackage::VPFloat& conj(const VPFloatPackage::VPFloat& x)  { return x; }
	inline const VPFloatPackage::VPFloat& real(const VPFloatPackage::VPFloat& x)  { return x; }
	inline VPFloatPackage::VPFloat imag(const VPFloatPackage::VPFloat&)    { return VPFloat(0.); }
	inline VPFloatPackage::VPFloat abs(const VPFloatPackage::VPFloat&  x)  { return abs(x); }
	inline VPFloatPackage::VPFloat abs2(const VPFloatPackage::VPFloat& x)  { return x*x; }
	inline VPFloatPackage::VPFloat sqrt(const VPFloatPackage::VPFloat& x)  { return VMath::vsqrt(x); }
}

inline bool operator== (const VPFloat& a, const int& b) { return a == VPFloat(b); }

namespace Eigen 
{

namespace internal 
{

inline bool isMuchSmallerThan(const VPFloatPackage::VPFloat& a, const VPFloatPackage::VPFloat& b, const VPFloatPackage::VPFloat& eps)
{
	return VPFloatPackage::abs(a) <= VPFloatPackage::abs(b) * eps;
}

inline bool isApprox(const VPFloatPackage::VPFloat& a, const VPFloatPackage::VPFloat& b, const VPFloatPackage::VPFloat& eps)
{
	return VPFloatPackage::VMath::isAlmostEqual(a,b,eps);
}

inline bool isApproxOrLessThan(const VPFloatPackage::VPFloat& a, const VPFloatPackage::VPFloat& b, const VPFloatPackage::VPFloat& eps)
{
	return a <= b or VPFloatPackage::VMath::isAlmostEqual(a,b,eps);
}

template<> inline double cast<VPFloatPackage::VPFloat,double>(const VPFloatPackage::VPFloat& x)
{ return x.getdouble(); }

template<> inline float cast<VPFloatPackage::VPFloat,float>(const VPFloatPackage::VPFloat& x)
{ return x.getfloat(); }

template<> inline long cast<VPFloatPackage::VPFloat,long>(const VPFloatPackage::VPFloat& x)
{ return long(x.getdouble()); }

template<> inline int cast<VPFloatPackage::VPFloat,int>(const VPFloatPackage::VPFloat& x)
{ return int(x.getdouble()); }

} // end namespace internal

}

#endif // EIGEN_VPFLOAT_SUPPORT_HPP
