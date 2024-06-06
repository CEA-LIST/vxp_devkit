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
 * Authors       : Alexandre Hoffmann, Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#include "bicgstab_kernel.hpp"
#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VBLAS.hpp"
#include "VPSDK/VMath.hpp" // for vsqrt

using namespace VPFloatPackage;

#define ITER_MAX 5


int bicgstab_vp(          // Solves Ax = b where A is a non-symetric square matrix using the Quasi-Minimal Residual method  
	int precision,          // precision used by VPFloat (aka the size of the mantissa)
	int transpose,
	int n,                  // size of the matrix 
	VPFloatArray& x,        // an initial guess for the system Ax = b
	matrix_t A,             // A
	matrix_t At,            // Transpose of A
	VPFloatArray b,         // LHS
	double tolerance,       // stoping criterion for the method. The method will stop if |r_k|^2 < tolerance^2 
	uint16_t exponent_size, // size of the exponent of a VPFloat
	int32_t stride_size)	    // mysterious variable which shall be 1
{
	short myBis=precision+exponent_size+1;
  int nbiter;

  VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
  VPFloatComputingEnvironment::set_precision(precision);
  VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);

  VPFloatArray x_k(exponent_size, myBis, stride_size, n );
  VPFloatArray r_k(exponent_size, myBis, stride_size, n );
  VPFloatArray hat_r_0(exponent_size, myBis, stride_size, n );
  
  VPFloatArray v_k(exponent_size, myBis, stride_size, n );
  VPFloatArray p_k(exponent_size, myBis, stride_size, n );
  VPFloatArray s_k(exponent_size, myBis, stride_size, n );
  VPFloatArray t_k(exponent_size, myBis, stride_size, n );
  
  VPFloatArray tmp(exponent_size, myBis, stride_size, n );
  VPFloatArray tmp2(exponent_size, myBis, stride_size, n );
  
  VPFloat squaredNorm_rk(exponent_size, myBis, stride_size );
  VPFloat r0_dot_vk(exponent_size, myBis, stride_size );
  VPFloat tk_dot_sk(exponent_size, myBis, stride_size );
  VPFloat squaredNorm_tk(exponent_size, myBis, stride_size );
  
  VPFloat rho(exponent_size, myBis, stride_size );
  VPFloat alpha(exponent_size, myBis, stride_size );
  VPFloat omega(exponent_size, myBis, stride_size );
  
  VPFloat rhoOld(exponent_size, myBis, stride_size );
  VPFloat beta(exponent_size, myBis, stride_size );
  
  /* x_k = {0} (choice) */
  VBLAS::vzero(precision, n, x_k);
  /* r_0 <- b - Ax_0 = b */
  VBLAS::vcopy(n, b, r_k);
  /* hat_r_0 <- r_0 (choice) */
  VBLAS::vcopy(n, r_k, hat_r_0);
  
  rho = 1.;
  alpha = 1.;
  omega = 1.;

  for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
		
		VBLAS::vdot(precision, n, r_k,  r_k, squaredNorm_rk);
		std::cout << "residus : " <<  double(squaredNorm_rk) << " - iter : " << nbiter << std::endl;
		if ((double)squaredNorm_rk < (tolerance*tolerance)) { VBLAS::vcopy(n, x_k,  x); break; }
		
		rhoOld = rho;
		VBLAS::vdot(precision, n, hat_r_0, r_k, rho);
		
		beta = (rho / rhoOld)*(alpha / omega);
		
		/* p <- r + beta*(p - omega*v) */
		VBLAS::vcopy(n, p_k, tmp);                    // tmp   = p
		VBLAS::vcopy(n, r_k, tmp2);                   // tmp2  = r
		VBLAS::vaxpy(precision, n, -omega, v_k, tmp); // tmp  -= omega*v  (tmp = p - omega*v)
		VBLAS::vaxpy(precision, n,  beta, tmp, tmp2); // tmp2 += beta*tmp (tmp = r + beta*(p - omega*v))
		VBLAS::vcopy(n, tmp2, p_k);
		/* v <- A*p */
		VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, 1.0, A, p_k, 0.0, v_k);
		VBLAS::vdot(precision, n, hat_r_0, v_k, r0_dot_vk);
		alpha = rho / r0_dot_vk;
		/* s <- r - alpha*v */
		VBLAS::vcopy(n, r_k, s_k);
		VBLAS::vaxpy(precision, n, -alpha, v_k, s_k);
		
		VBLAS::vdot(precision, n, s_k,  s_k, squaredNorm_rk);
		if ((double)squaredNorm_rk < (tolerance*tolerance)) { VBLAS::vcopy(n, x_k,  x); break; }

		/* t <- A*s */
		VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, 1.0, A, s_k, 0.0, t_k);
		/* omega <- (t,s) / (t,t) */
		VBLAS::vdot(precision, n, t_k, s_k, tk_dot_sk);
		VBLAS::vdot(precision, n, t_k, t_k, squaredNorm_tk);
		omega = tk_dot_sk / squaredNorm_tk;
		/* x <- x + alpha*p + omega*s */
		VBLAS::vaxpy(precision, n,  alpha, p_k, x_k);
		VBLAS::vaxpy(precision, n,  omega, s_k, x_k);
		/* r <- s - omega*t */
		VBLAS::vcopy(n, s_k, r_k);                    // r  = s
		VBLAS::vaxpy(precision, n, -omega, t_k, r_k); // r -= omega*t
	}
  if (nbiter==(n*ITER_MAX)) // ca n'a pas converge
  {
		nbiter=-2; //VERYGRANDNOMBRE-1;
	}
	
  return (nbiter + 1);
}
