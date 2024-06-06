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

#include "qmr_kernel.hpp"
#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VBLAS.hpp"
#include "VPSDK/VMath.hpp" // for vsqrt

using namespace VPFloatPackage;

#define ITER_MAX 5

/*
 * ------------------- QMR ---------------------
 * Freund, R. W., & Nachtigal, N. M. (1994). An implementation of the QMR method based on coupled two-term recurrences. SIAM Journal on Scientific Computing, 15(2), 313-337. (Algorithm 7.3)
 * Barrett, R., Berry, M., Chan, T. F., Demmel, J., Donato, J., Dongarra, J., ... & Van der Vorst, H. (1994). Templates for the solution of linear systems: building blocks for iterative methods. Society for Industrial and Applied Mathematics.
 */

int qmr_vp(               // Solves Ax = b where A is a non-symetric square matrix using the Quasi-minimal residual method from [1,2]    
	int precision,          // precision used by VPFloat (aka the size of the mantissa)
	int transpose,			// flag used to specify that matrices are transposed
	int n,                  // size of the matrix 
	VPFloatArray& x,        // an initial guess for the system Ax = b
	matrix_t A,             // A
	matrix_t At,            // Transpose of A
	VPFloatArray b,         // LHS
	double tolerance,       // stoping criterion for the method. The method will stop if |r_k|^2 < tolerance^2 
	uint16_t exponent_size, // size of the exponent of a VPFloat
	int32_t stride_size)    // mysterious variable which shall be 1
{
	//inline void norm(int precision, int n, VPFloatArray vec, VPFloat& res) { VBLAS::vdot(precision, n, vec,  vec, res); res = VMath::vsqrt(res); }
	
	short myBis=precision+exponent_size+1;

  VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
  VPFloatComputingEnvironment::set_precision(precision);
  VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);

  int nbiter;
  VPFloatArray x_k(exponent_size, myBis, stride_size, n );
  VPFloatArray r_k(exponent_size, myBis, stride_size, n );
  
  VPFloatArray tilde_v_k(exponent_size, myBis, stride_size, n ); // un-normalized version of v_k
  VPFloatArray tilde_w_k(exponent_size, myBis, stride_size, n ); // un-normalized version of w_k
  
  VPFloatArray v_k(exponent_size, myBis, stride_size, n );
  VPFloatArray w_k(exponent_size, myBis, stride_size, n );
  
  VPFloatArray p_k(exponent_size, myBis, stride_size, n );
  VPFloatArray q_k(exponent_size, myBis, stride_size, n );
  VPFloatArray d_k(exponent_size, myBis, stride_size, n );
  VPFloatArray s_k(exponent_size, myBis, stride_size, n );
  
  VPFloatArray Ap_k(exponent_size, myBis, stride_size, n );
  VPFloatArray Atq_k(exponent_size, myBis, stride_size, n );
  
  VPFloat beta_k(exponent_size, myBis, stride_size );
  VPFloat gamma_k(exponent_size, myBis, stride_size );
  
  VPFloat c_k(exponent_size, myBis, stride_size );
  VPFloat mu_k(exponent_size, myBis, stride_size );
  VPFloat lambda_k(exponent_size, myBis, stride_size );
  VPFloat sigma_k(exponent_size, myBis, stride_size );
  VPFloat varTheta_k(exponent_size, myBis, stride_size );
  VPFloat eta_k(exponent_size, myBis, stride_size );
  
  VPFloat beta_kp1(exponent_size, myBis, stride_size );
  VPFloat gamma_kp1(exponent_size, myBis, stride_size );
  
  VPFloat c_km1(exponent_size, myBis, stride_size );
  VPFloat varTheta_km1(exponent_size, myBis, stride_size );
  
  VPFloat ONE(exponent_size, myBis, stride_size ); ONE = 1;
  VPFloat squaredNorm_rk(exponent_size, myBis, stride_size );
  
  VPFloat squaredNorm_q(exponent_size, myBis, stride_size );

  /* x_k = {0} (choice) */
  VBLAS::vzero(precision, n, x_k);
  /* r_0 <- b - Ax_0 = b */
  /* tilde_v_0 <- r_0 */
  /* tilde_w_0 <- r_0 */
  VBLAS::vcopy(n, b, r_k);
  VBLAS::vcopy(n, r_k, tilde_v_k);
  VBLAS::vcopy(n, r_k, tilde_w_k);
  /* beta_1  <- |tilde_v_0| */
  /* gamma_1 <- |tilde_w_0| = |tilde_v_0| = beta_1 */
  VBLAS::vdot(precision, n, tilde_v_k, tilde_v_k, beta_k); beta_k = VMath::vsqrt(beta_k);
  gamma_k = beta_k;
  /* p_0 <- q_0 <- d_0 <- s_0 <- {0} for the first iteration*/
  VBLAS::vzero(precision, n, p_k);
  VBLAS::vzero(precision, n, q_k);
  VBLAS::vzero(precision, n, d_k);
  VBLAS::vzero(precision, n, s_k);
  /* c_0 <- mu_0 <- 1 */
  /* varTheta_0 <- 0 */
  /* eta_0 <- -1 */
  c_k = 1;
  mu_k = 1;
  varTheta_k = 0;
  eta_k = -1;
  
  for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
		c_km1 = c_k;
		varTheta_km1 = varTheta_k;

		VBLAS::vdot(precision, n, r_k,  r_k, squaredNorm_rk);
		std::cout << "residus : " <<  double(squaredNorm_rk) << " - iter : " << nbiter << std::endl;
		if ((double)squaredNorm_rk < (tolerance*tolerance)) { VBLAS::vcopy(n, x_k,  x); break; }
		/* =========================== */
		/* Coupled two term recurrence */
		/* =========================== */
		VBLAS::vcopy(n, tilde_v_k, v_k);
		VBLAS::vcopy(n, tilde_w_k, w_k);
		VBLAS::vscal(precision, n, ONE / beta_k, v_k);
		VBLAS::vscal(precision, n, ONE / gamma_k, w_k);
		/* sigma_k <- (w_k, v_k) */
		VBLAS::vdot(precision, n, w_k, v_k, sigma_k);

		/* p_k <- v_k - p_{k-1}*(gamma_k*sigma_k / mu_{k-1}) */
		VBLAS::vscal(precision, n, -(gamma_k*sigma_k / mu_k), p_k); // here mu_k hasn't been updated yet, mu_k thus contains the value of mu_{k-1}
		VBLAS::vaxpy(precision, n, ONE, v_k, p_k);
		/* q_k <- w_k - q_{k-1}*(beta_k*sigma_k / mu_{k-1})  */
		VBLAS::vscal(precision, n, -(beta_k*sigma_k / mu_k), q_k); // here mu_k hasn't been updated yet, mu_k thus contains the value of mu_{k-1}
		VBLAS::vaxpy(precision, n, ONE, w_k, q_k);
		/* Compute Ap_k and At q_k */
		VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, ONE, A, p_k, 0.0, Ap_k);
		VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, ONE, At, q_k, 0.0, Atq_k);
		/* mu_k <- (q_k, Ap_k) */
		VBLAS::vdot(precision, n, q_k, Ap_k, mu_k);
		/* lambda_k <- mu_k / sigma_k */
		lambda_k = mu_k / sigma_k;

		/* tilde_v_{k+1} <- Ap_k - lambda_k v_k */
		VBLAS::vcopy(n, Ap_k, tilde_v_k);
		VBLAS::vaxpy(precision, n, -lambda_k, v_k, tilde_v_k);
		/* tilde_w_{k+1} <- At q_k - lambda_k w_k */
		VBLAS::vcopy(n, Atq_k, tilde_w_k);
		VBLAS::vaxpy(precision, n, -lambda_k, w_k, tilde_w_k);
		/* beta_{k+1} <- \| tilde_v_{k+1} \| */
		/* gamma_{k+1} <- \| tilde_w_{k+1} \| */
		VBLAS::vdot(precision, n, tilde_v_k, tilde_v_k, beta_kp1); beta_kp1 = VMath::vsqrt(beta_kp1);
		VBLAS::vdot(precision, n, tilde_w_k, tilde_w_k, gamma_kp1); gamma_kp1 = VMath::vsqrt(gamma_kp1);

		/* ====================== */
		/* Quasi minimal residual */
		/* ====================== */
		varTheta_k = beta_kp1 / (c_km1*VPFloatPackage::abs(lambda_k));

		c_k = ONE / VMath::vsqrt(ONE + varTheta_k*varTheta_k);
		eta_k *= -((beta_k*c_k*c_k) / (lambda_k*c_km1*c_km1));
		/* d_k <- eta_k p_k + (varTheta_{k-1}*c_k)^2 d_{k-1} */
		VBLAS::vscal(precision, n, (varTheta_km1*c_k)*(varTheta_km1*c_k), d_k);
		VBLAS::vaxpy(precision, n, eta_k, p_k, d_k);
		/* s_k <- eta_k Ap_k + (varTheta_{k-1}*c_k)^2 s_{k-1} */
		VBLAS::vscal(precision, n, (varTheta_km1*c_k)*(varTheta_km1*c_k), s_k);
		VBLAS::vaxpy(precision, n, eta_k, Ap_k, s_k);
		/* x_{k} <- x_{k-1} + d_k */
		/* r_{k} <- r_{k-1} - s_k */
		VBLAS::vaxpy(precision, n,  ONE, d_k, x_k);
		VBLAS::vaxpy(precision, n, -ONE, s_k, r_k);

		/* ======================== */
		/* Setup for next iteration */
		/* ======================== */
		beta_k = beta_kp1;
		gamma_k = gamma_kp1;
	}
  if (nbiter==(n*ITER_MAX)) // ca n'a pas converge
  {
		nbiter=-2; //VERYGRANDNOMBRE-1;
	}
	
  return (nbiter + 1);
}
