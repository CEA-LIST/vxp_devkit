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
 * Authors       : Yves Durand, Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

/* 
 * gestion de la structure crs 
 * == matrices creuses
 */
#include "cg_kernel.hpp"
#include "VPSDK/VBLAS.hpp"
#include "VRPSDK/perfcounters/cpu.h"

// package de support VPFloat
using namespace VPFloatPackage;

/* 
 * les paramètres
 * ==========================
 */
// ITER_MAX: on abandonne après ndiag*ITER_MAR iterations
//#define ITER_MAX 50
#define ITER_MAX 5


/*
 * ----------------------------------
 * ----dense CG ---------------------
 */
int cg_vp(int precision,
    int transpose,
    int n,
    VPFloatArray & x,  // valeur de sortie 
    matrix_t A,
    VPFloatArray & b,
    double tolerance,
    uint16_t exponent_size,
    int32_t stride_size)
{

  short myBis=precision+exponent_size+1;
  int nbiter;


  VPFloatComputingEnvironment::set_precision(myBis);
  VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);

  VPFloatArray r_k(exponent_size, myBis, stride_size, n );
  VPFloatArray p_k(exponent_size, myBis, stride_size, n );
  VPFloatArray Ap_k(exponent_size, myBis, stride_size, n );
  VPFloatArray x_k(exponent_size, myBis, stride_size, n );
  VPFloat      alpha   (exponent_size, myBis, stride_size );
  VPFloat      beta    (exponent_size, myBis, stride_size );
  VPFloat      rs      (exponent_size, myBis, stride_size );
  VPFloat      rs_next (exponent_size, myBis, stride_size );

  /* r_k <- b */
  VBLAS::vcopy(n,b,r_k);
#ifdef DBG
  v_disp_matrix(r_k,"r_k",n,1);
#endif
  /* p_k <- b */
  VBLAS::vcopy(n,b,p_k);
  /* x_k = {0} */
  VBLAS::vzero(precision, n, x_k);
  VBLAS::vzero(precision, n, Ap_k);
  rs_next = 0.0;
 
  // rs = rk'*rk
  VBLAS::vdot(precision, n, r_k, r_k, rs);
  
  for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
    std::cout << "residus : " <<  double(rs_next) << " - iter : " << nbiter << std::endl;

    // Ap_k = A * p_k
    VBLAS::vgemvd(precision,
    transpose == 0 ? 'N' : 'Y',
		n, //m
		n, //n
		1.0, //CONST_1, //alpha 
		A,
		p_k /* x */,
		0.0, //CONST_0 /*beta*/,
		Ap_k  /*Y*/);
#ifdef DBG
    v_disp_matrix(Ap_k,"Ap_k",n,1);
#endif
    
    // alpha = Ap_k * p_k
    VBLAS::vdot(precision, n, p_k, Ap_k,  alpha);
#ifdef DBG
    std::cout << "iter "<<nbiter<<" |p_k A p_k|="<<alpha.getdouble()<<"\n";
#endif

    // alpha = rs/alpha
    alpha = rs/alpha;
#ifdef DBG
    std::cout << "iter "<<nbiter<<" \nalpha="<<alpha.getdouble()<<"\n";
#endif

    // x_k = x_k + alpha*p_k
    VBLAS::vaxpy(precision, n, alpha, p_k,  x_k);

#ifdef DBG
    v_disp_matrix(x_k,"x_k = x_k+alpha*p_k",n,1);
#endif
	
    // r_k = r_k - alpha * Ap_k

#ifdef DBG
    std::cout<<"alpha = "<< alpha <<" \n";
#endif

    VBLAS::vaxpy(precision, n, -alpha, Ap_k, r_k);
#ifdef DBG
    v_disp_matrix(r_k,"r_k = r_k-alpha*Ap_k",n,1);
#endif

    // rs_next (rs:r square)
    //double rs_sqrt;
    VBLAS::vdot(precision, n, r_k,  r_k, rs_next);
    
#ifdef DBG
    {
      printf("\t(debug) iter : %d - rs_next = %f\n", nbiter, (double)rs_next);
      std::cout <<" rs_next :" <<rs_next<<"\n";
    }
#endif
#if DEBUG_DATA_CG
      {
	vprint_vector("\t(debug) x_k", n, x_k, ENV);
      }
#endif
      if ((double)rs_next < (tolerance*tolerance)) {
	VBLAS::vcopy(n, x_k,  x);
	break;
      }
      // p_k = p_k * (rs_next/rs)
      // reutilisons rs pour beta
      beta= rs_next/rs;
      
      VBLAS::vscal(precision, n, beta, p_k);

      // p_k = p_k + r_k
      VBLAS::vaxpy(precision, n, 1.0 /*CONST_1*/, r_k,  p_k);

      // rs = rs_next
      rs = rs_next;
    }
  if (nbiter==(n*ITER_MAX))
    // ca n'a pas converge
    nbiter=-2; //VERYGRANDNOMBRE-1;

  return (nbiter + 1);
}

