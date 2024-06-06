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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "MTXUtil/crs.h"

//#define _MVM(M,V,R) fastCrsMVM((M),(V),(R))

//
/*
 *** acces au elements de la matrice A[row,col]
 * attention,  optimisé pour un ordre de rang croissants
 * A UTILISER PARCIMONIEUSEMENT
 */

double crsget(crs_t * spm , int row, int col)  
{
  double val =0.0;
  // debut de la ligne. Tention, row est 1:N
  int firstOfRow = (spm->ROW_PTR[row-1])-1;
  int nextRow    = (spm->ROW_PTR[row])-1;
  int idx;
  for (idx=firstOfRow; idx<nextRow; idx++) {
    if (spm->COL_IND[idx]==col) {
      val = spm->VAL[idx];
      break;
    }
  } 
  return(val);
}

double crsgetC(crs_t * spm , int row, int col)  
{
  return(crsget(spm,row+1,col+1));
}

void prettyCrs(char * name, crs_t *spm) {
  // inefficace, test de crsget
  int i,j;
  printf("%s=[",name);
  for (i=1; i<=spm->M; i++) {
    if (i>1) printf(";\n");
    for (j=1; j<=spm->N; j++) {
      if (i==spm->M && j==spm->N) {
	printf(" %f ]\n",crsget(spm,i,j));
      } else { printf(" %f ",crsget(spm,i,j)); }
    }
  }
}

int CrsIsSymmetric(crs_t * spm){
  int i,j;

  fprintf(stderr,"checking for Symmetry \n");
  if (spm->M!=spm->N) {
    fprintf(stderr,"SptIsSymmetric: Come on, Matrix is not square!\n");
    return(-1);
  }
  for (i=1; i<=spm->M; i++) {
    for (j=1; j<=spm->N; j++) {
      //fprintf(stderr,"..%d,%d..  ",i,j);

      if (crsget(spm,i,j) != crsget(spm,i,j)) {
	printf(" symmetry mismatch @ A[%d,%d] != A[%d,%d]\n",
	       i,j,j,i);
	return(-1);
      }
    }}
  return(1);
}


/* multiplication matrice creuse x vecteur dense
   NE PAS UTILISER pour le calcul
   uniquement POUR REFERENCE et non-regression
   la bonne methode c'est fastCrsMVM
*/
int slowCrsMVM (crs_t * A, double* V, double *RES) {
  // Res=A*V'
  // convention C: A[0:m-1; 0:n-1]
  // sanity check
  int n,m;
  n=A->N; m=A->M;
  int i,k;
  for (i=0; i<m; i++) {
    RES[i]=0.0;
    for (k=0; k<n; k++) {
      RES[i]=RES[i]+crsgetC(A,i,k) * V[k];
    }
  }
  return(1);
}

int fastCrsMVM (crs_t * A, double* x, double *y) {
  // y=A*x
  // algo assez evident mais analyse ds
  // "Cache oblivious sparse matrix-vector multiplication
  int m=A->M;
  int j, i,k;
  for (i=0; i<m; i++) {
    y[i]=0.0;
    for (k=(A->ROW_PTR[i])-1; k<(A->ROW_PTR[i+1])-1; k++) {
      // k est l'indice de ligne
      j = (A->COL_IND[k])-1; // tjrs le decalage magique
      y[i]=y[i]+(A->VAL[k])*x[j];
      //printf("\n k=%d  A->VAL[k]=%f  y[%d]=%f ,",k,A->VAL[k],i,y[i]);
    }
  }
  return(1);
  
}



void dumpCrs(crs_t * spm) {
  int i;
  printf("ROW_PTR :\n");
  for (i=0; i<spm->N; i++) printf("%d ",spm->ROW_PTR[i]);
  printf("\n VAL :\n");    
  for (i=0; i<spm->NZNUM; i++) printf("%f ",spm->VAL[i]);
  printf("\n COL_IND :\n");    
  for (i=0; i<spm->NZNUM; i++) printf("%d ",spm->COL_IND[i]);
  printf("\n");
}

void crs2dense(crs_t *spm, double * M){
  int i, j, ii;
  for (ii=0; ii<spm->M*spm->N; ii++) M[ii]=0.0;
  for (i=1; i<=spm->M; i++) {
    for (j=1; j<=spm->N; j++) {
      M[(i-1)*spm->N+(j-1)]=crsget(spm,i,j);
      //printf("M[%d,%d]= %f \n",i-1,j-1,crsget(spm,i,j));
    }
  }
}
int crsIdent(int n, crs_t * Ident) {
  int *myROW_PTR, *myCOL_IND;
  int ii;
  double *myval;
  myROW_PTR=(int *)malloc((n+1)*sizeof(int));
  myCOL_IND=(int *)malloc(n*sizeof(int));
  myval = (double *)malloc(n*sizeof(double));
  for (ii=0; ii<n; ii++) {
    myval[ii]=(double) 1.0; //ii+1.0; // debug!
    myROW_PTR[ii]=ii+1;
    myCOL_IND[ii]=ii+1;
  }
  myROW_PTR[n]=n+1;
  Ident->ROW_PTR=myROW_PTR;
  Ident->COL_IND=myCOL_IND;
  Ident->VAL=myval;
  Ident->NZNUM=n;
  Ident->M=n;
  Ident->N=n;
  return(n);
}

void crsFree(struct crs *spm) {
  if (spm->ROW_PTR  != NULL ) free(spm->ROW_PTR);
  if (spm->COL_IND  != NULL ) free(spm->COL_IND);
  if (spm->VAL  != NULL ) free(spm->VAL);
  free(spm);
}