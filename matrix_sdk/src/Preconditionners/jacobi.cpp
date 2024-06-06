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

#include <Preconditionners.hpp>
#include <OSKIHelper.hpp>
#include <iostream>
#include "MTXUtil/crs.h"

matrix_t jacobi(matrix_t a_input_matrix, double a_shifter) {

    if ( a_input_matrix->type_value == COMPLEX_VALUE ) {
        std::cout << __FUNCTION__ << "Not implemented for complex values." << std::endl;
        return NULL;
    }

    double * l_diag = (double *)VPFloatPackage::OSKIHelper::getDiag(a_input_matrix);

    // *** creation de la structure i,j,val a plat
    // *** semble inutile: cela l'est effectivement
    int nz=a_input_matrix->n;   
    int *ROW, *COL;
    double *val;
    ROW = (int *) malloc(nz * sizeof(int));
    COL = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    for (int ii=0; ii<a_input_matrix->n; ii++) {
        ROW[ii]=ii+1;
        COL[ii]=ii+1;
        if ( l_diag[ii] == 0.0 ) {
            std::cout << __FUNCTION__ << " Diagonal element " << ii << " of matrix is null!! Jacobi may contain NaN!" << std::endl;
            printf(" Diagonal element %d of matrix is null!! Jacobi may contain NaN!", ii);
        }
        val[ii]=1.0/(l_diag[ii] + a_shifter);
         
    }

    // *** creation de la struc crs SPM
    struct crs * SPM;
    SPM=(struct crs *)malloc(sizeof(struct crs));
    SPM->M=a_input_matrix->n;
    SPM->N=a_input_matrix->n;
    SPM->NZNUM=a_input_matrix->n;
    SPM->VAL=val;
    SPM->COL_IND =(int *) malloc(SPM->NZNUM*sizeof(int));
    SPM->ROW_PTR =(int *) malloc((SPM->M+1)*sizeof(int));

    int idx;       // parcourt la triple liste _VAL, _ROW, _COL
    int new_row=1, previous_row=0; // logic
    SPM->ROW_PTR[0]=ROW[0];
    for (idx=0; idx<SPM->NZNUM; idx++) {
        SPM->COL_IND[idx]=COL[idx];
        if (ROW[idx]==previous_row) {
            // skip
        } else {
            new_row=ROW[idx];
            // attention ! les valeurs de ROW_PTR sont de 1 à N -> add 1
            // le stockage ds ROW_PTR sont de 0 a M-1 -> retrancher 1 !!
            for (int kk=previous_row; kk <= (new_row-1); kk++) {
                SPM->ROW_PTR[kk]=idx+1;
            }
            previous_row=new_row;
        }
        SPM->ROW_PTR[a_input_matrix->n]=a_input_matrix->n+1; // convention IMPORTANTE
    }
    matrix_t l_iM=buildCSR(a_input_matrix->n,a_input_matrix->n,SPM->ROW_PTR,SPM->COL_IND,SPM->VAL,1); // 1==ONE_BASED

    return l_iM;
}