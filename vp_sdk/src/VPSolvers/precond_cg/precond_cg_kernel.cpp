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
#include "precond_cg_kernel.hpp"
#include "VPSDK/VBLAS.hpp"
#include <cmath>

// package de support VPFloat
using namespace VPFloatPackage;

#define ITER_MAX 10


/* ----------------------------------
 * ---- sparse precond CG -----------
 * ----------------------------------  */
int precond_cg_vp(  int precision,
                    int transpose,
                    int n,
                    VPFloatArray & x,  // valeur de sortie 
                    matrix_t A,
                    matrix_t iM,
                    VPFloatArray b,
                    double tolerance,
                    uint16_t exponent_size,
                    int32_t stride_size)

/*
 * inspired by Saad "Iterative Methods ..." algorithm 9.1
 * en fixant x0=[0 ... 0]
 */
{
    short myBis=precision+exponent_size+1;

    VPFloatComputingEnvironment::set_rounding_mode(VP_RNE);
    VPFloatComputingEnvironment::set_precision(precision);
    VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);
//     VPFloatComputingEnvironment::set_tempory_var_environment(VPFLOAT_EVP_MAX.es, VPFLOAT_EVP_MAX.bis, VPFLOAT_EVP_MAX.stride);

    int nbiter;
    VPFloatArray r_j(exponent_size, myBis, stride_size, n );
    VPFloatArray z_j(exponent_size, myBis, stride_size, n );
    VPFloatArray p_j(exponent_size, myBis, stride_size, n );
    VPFloatArray Ap_j(exponent_size, myBis, stride_size, n );
    VPFloatArray x_j(exponent_size, myBis, stride_size, n );
    VPFloat      alpha_j (exponent_size, myBis, stride_size );
    VPFloat      beta_j  (exponent_size, myBis, stride_size );
    VPFloat      r_jxz_j (exponent_size, myBis, stride_size );
    VPFloat      r_jxz_jnext (exponent_size, myBis, stride_size );
    VPFloat      r_jsq  (exponent_size, myBis, stride_size );
    VPFloat      Apjxpj (exponent_size, myBis, stride_size );
    VPFloat      rs_next (exponent_size, myBis, stride_size );

    /* r_j <- b */
    VBLAS::vcopy(n,b,r_j);

#ifdef DBG
    v_disp_matrix(r_j,"r_j",n,1);
#endif
    /* z_j=iM * r_j */
    VBLAS::vgemvd(  precision, transpose == 0 ? 'N' : 'Y',
                    n,   //m
                    n,   //n
                    1.0, //alpha 
                    iM,
                    r_j /* x */,
                    0.0, //CONST_0 /*beta*/,
                    z_j  /*Y*/);

    /* p_j <- z_j */
    VBLAS::vcopy(n,z_j,p_j);

    /* x_k = {0} */
    VBLAS::vzero(precision, n, x_j);
    VBLAS::vzero(precision, n, Ap_j);
    //  rj_next = 0.0;
    VBLAS::vdot(precision, n, r_j, z_j,  r_jxz_j);

    for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
        std::cout << "residus : " <<  double(r_jsq) << " - iter : " << nbiter << std::endl;

        // Ap_j = A * p_j
        VBLAS::vgemvd(  precision, transpose == 0 ? 'N' : 'Y',
                        n,   //m
                        n,   //n
                        1.0, //alpha 
                        A,
                        p_j /* x */,
                        0.0, //CONST_0 /*beta*/,
                        Ap_j  /*Y*/);

#ifdef DBG
        v_disp_matrix(Ap_j,"Ap_j",n,1);
#endif

        // \alpha_j = (r-j,z_j) / (Ap_j,p_j)
        // d'abord Apjxpj = Ap_j * p_j
        VBLAS::vdot(precision, n, p_j, Ap_j,  Apjxpj);

#ifdef DBG
        std::cout << "iter "<<nbiter<<" |p_j A p_j|="<<(double)alpha<<"\n";
#endif
        // ensuite r_j x z_j
        // ATTENTION  A OPTIMISER
        //VBLAS::vdot(precision, n, r_j, z_j,  r_jxz_j);
      
        // alpha = r_jxz_j/Apjxpj
        
        // VPFloatComputingEnvironment::set_tempory_var_environment(VPFLOAT_EVP_MAX.es, VPFLOAT_EVP_MAX.bis, VPFLOAT_EVP_MAX.stride);
        // printf("r_jxz_j: %e - Apjxpj:%e\n", double(r_jxz_j), double(Apjxpj));
        alpha_j = r_jxz_j/Apjxpj;
        // printf("alpha_j: %e\n", double(alpha_j));

        // VPFloatComputingEnvironment::set_tempory_var_environment(exponent_size, myBis, 1);
#ifdef DBG
        std::cout << "iter "<<nbiter<<" \nalpha="<<alpha<<"\n";
#endif

        // x_j = x_j + alpha*p_j
        VBLAS::vaxpy(precision, n, alpha_j, p_j,  x_j); 
        // for (int ii=0; ii<n; ii++) x_j[ii]=x_j[ii]+alpha_j*p_j[ii];
#ifdef DBG
        v_disp_matrix(x_k,"x_k = x_k+alpha_j*p_j",n,1);
#endif

        // r_{j+1} = r_j - alpha_j * Ap_j
        // la variable r_j devient r_{j+1}

#ifdef DBG
        Ap_j.printAsVector("Ap_j:");
#endif
        
        VBLAS::vaxpy(precision, n, -alpha_j, Ap_j, r_j);
#ifdef DBG
        v_disp_matrix(r_j,"r_j = r_j-alpha*Ap_j",n,1);
#endif
        // calcul de r_jsq
        VBLAS::vdot(precision, n, r_j, r_j,  r_jsq);

        // printf("r_jsq: %e\n", double(r_jsq));

        //if (vtod(rs_next, ENV) < (tolerance*tolerance)) {
        if (r_jsq < (tolerance*tolerance)) {
            VBLAS::vcopy(n, x_j,  x);
            break;
        }
        // if ( ( double(alpha_j) == INFINITY)  || std::isnan(double(r_jsq)) ) {
        //         printf("Test to remove reached.\n");
        //         break;
        // }
        // z_{j+1}=iM * r_{j+1}
        // on remplace z_j
        VBLAS::vgemvd(  precision, transpose == 0 ? 'N' : 'Y',
                        n,   //m
                        n,   //n
                        1.0, //alpha 
                        iM,
                        r_j /* x */,
                        0.0, //CONST_0 /*beta*/,
                        z_j  /*Y*/);
        // r_jxz_jnzext
        // r_jxz_j a encore la val precedente
        VBLAS::vdot(precision, n, r_j, z_j,  r_jxz_jnext);

        beta_j= r_jxz_jnext/r_jxz_j;

        VBLAS::vscal(precision, n, beta_j, p_j);

        // p_j = z_j + r_j
        VBLAS::vaxpy(precision, n, 1.0 /*CONST_1*/, z_j,  p_j);

        // rs = rs_next
        r_jxz_j = r_jxz_jnext;

        
        // printf("nbiter:%d\n", nbiter);
        
    }

    if (nbiter==(n*ITER_MAX)) {
        // ca n'a pas converge
        nbiter=-2; //VERYGRANDNOMBRE-1;
    }

    return (nbiter + 1);
}
