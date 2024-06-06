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

#include "bicg_kernel.hpp"
#include "VPSDK/VBLAS.hpp"

using namespace VPFloatPackage;
/* 
 * les paramètres
 * ==========================
 */
// ITER_MAX: on abandonne après ndiag*ITER_MAR iterations
//#define ITER_MAX 50
#define ITER_MAX 5
//#define DBG

/*
 * ------------------- BiCG ---------------------
 * Saad algorithme 7.3
 */
int bicg_vp(int precision,
      int transpose,
	    int n,
	    VPFloatArray & x,  // valeur de sortie et d'entree
	    matrix_t A,
	    matrix_t At,  // petite arnaque en attendant le support vgemv
	    VPFloatArray b,
	    double tolerance,
      uint16_t exponent_size,
      int32_t stride_size)
{
  short myBis=precision+exponent_size+1;
  int nbiter;

  VPFloatComputingEnvironment::set_precision(precision);
  VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);

  VPFloatArray r_k(exponent_size, myBis, stride_size, n );
  VPFloatArray rstar_k(exponent_size, myBis, stride_size, n );
  VPFloatArray p_k(exponent_size, myBis, stride_size, n );
  VPFloatArray pstar_k(exponent_size, myBis, stride_size, n );
  VPFloatArray Ap_k(exponent_size, myBis, stride_size, n );
  VPFloatArray minus_alphaxAtpstar_k
                   (exponent_size, myBis, stride_size, n );
  VPFloatArray x_k (exponent_size, myBis, stride_size, n );
  VPFloat      alpha   (exponent_size, myBis, stride_size );
  VPFloat      alpha_denom (exponent_size, myBis, stride_size );
  VPFloat      beta    (exponent_size, myBis, stride_size );
  VPFloat      rs      (exponent_size, myBis, stride_size );
  VPFloat      rxrstar (exponent_size, myBis, stride_size );
  VPFloat      rxrstar_next
                       (exponent_size, myBis, stride_size );
  VPFloat      rs_next (exponent_size, myBis, stride_size );;

  VPFloatArray Ax_check(exponent_size, myBis, stride_size, n );
  VPFloatArray Ax_minusB_check (exponent_size, myBis, stride_size, n );
  VPFloat      norm_Ax_minusB_check (exponent_size, myBis, stride_size );

  /* x_k = {0} (choix) */
  VBLAS::vzero(precision, n, x_k);
  /* r_k <- b - Ax0 */
  /* on considere x0 ca fait r_k <- b */
  VBLAS::vcopy(n,b,r_k);
  /* r*_k <- r0 (choix, ok classique selon Saad) */
  VBLAS::vcopy(n,b,rstar_k);
  /* p_k <- b */
  VBLAS::vcopy(n,r_k    ,p_k);
  VBLAS::vcopy(n,rstar_k,pstar_k);
  VBLAS::vcopy(n,b,Ax_minusB_check);

  rs_next = 0.0;
  // rxrstar = rk'*rstark
  VBLAS::vdot(precision, n, r_k, rstar_k, rxrstar);

  for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
    //std::cout<<"------ITER "<<nbiter<<" ----------\n";
		std::cout << "residus : " <<  double(rs_next) << " - iter : " << nbiter << std::endl;
    if ( rs_next.isNaN() ) {
      std::cout << "rs_next is NaN. Abort at solver at iteration " << nbiter << " ! " << std::endl;
      nbiter = -2;
      break;
    }
    // Ap_k = A * p_k
    VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, 1.0, //CONST_1, //alpha 
		  A, p_k, // x
		  0.0, //CONST_0 /*beta*/,
		  Ap_k  /*Y*/);

    // alphadenom = pstar_k' * Ap_k
    VBLAS::vdot(precision, n, pstar_k, Ap_k,  alpha_denom);

    // alpha = rs/alphadenom
    alpha = rxrstar/alpha_denom;
#ifdef DBG
    std::cout << "iter "<<nbiter<<"alpha="<<(double)alpha<<"\n";
#endif

    // x_k = x_k + alpha*p_k
    VBLAS::vaxpy(precision, n, alpha, p_k,  x_k);
	
    // r_k = r_k - alpha * Ap_k
    VBLAS::vaxpy(precision, n, -alpha, Ap_k, r_k);

    // rs_next (rs:r square)
    //double rs_sqrt;
    VBLAS::vdot(precision, n, r_k,  r_k, rs_next);
    
#ifdef DBG
    {
      printf("\t(debug) rs_next = %f\n", (double)rs_next);
    }
#endif
#ifdef DBG
      {
	vprint_vector("\t(debug) x_k", n, x_k, ENV);
      }
#endif



      if ((double)rs_next < (tolerance*tolerance)) {
	VBLAS::vcopy(n, x_k,  x);
	break;
      }
      // rstar_k=rstar_k - alpha_k A^T pstar_k
      // en 2 etapes
      VBLAS::vgemvd(precision, transpose == 0 ? 'N' : 'Y', n, n, 1.0, //-alpha, 
		    At, pstar_k, // x
		    0.0, //CONST_0 /*beta*/,
		    minus_alphaxAtpstar_k  /*Y*/);
      VBLAS::vscal(precision, n, -alpha, minus_alphaxAtpstar_k);

      VBLAS::vaxpy(precision, n, 1.0 /*CONST_1*/, minus_alphaxAtpstar_k,  rstar_k);
      // reutilisons rs pour beta
      // rxrstar_next = rk'*rstark
      VBLAS::vdot(precision, n, r_k, rstar_k, rxrstar_next);

      beta= rxrstar_next/rxrstar;
      
      VBLAS::vscal(precision, n, beta, p_k);
      // p_k = p_k + r_k
      VBLAS::vaxpy(precision, n, 1.0 /*CONST_1*/, r_k,  p_k);

      VBLAS::vscal(precision, n, beta, pstar_k);

      // pstar_k = p_k + r_k
      VBLAS::vaxpy(precision, n, 1.0 /*CONST_1*/, rstar_k,  pstar_k);

      // rs = rs_next
      rxrstar = rxrstar_next;

    }
  if (nbiter==(n*ITER_MAX)) 
    // ca n'a pas converge
    nbiter=-2; //VERYGRANDNOMBRE-1;
    
  if ( nbiter > 0 ) {
    // Compute ||Ax-b||
    VBLAS::vgemvd(precision,
                  transpose == 0 ? 'N' : 'Y',
                  n,
                  n,
                  1.0, //CONST_1, //alpha 
		              A,
                  x, // x
		              0.0, //CONST_0 /*beta*/,
		              Ax_check  /*Y*/);

    VBLAS::vaxpy(precision, n, -1.0, Ax_check,  Ax_minusB_check);

    VBLAS::vdot(precision, n, Ax_minusB_check, Ax_minusB_check, norm_Ax_minusB_check);

    std::cout << "||Ax-b|| : " << (double)norm_Ax_minusB_check << " - tolerance : " << tolerance << std::endl;

  }

  return (nbiter + 1);
}