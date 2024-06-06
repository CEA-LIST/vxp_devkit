
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

#ifndef __YDCRS__
#define __YDCRS__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { REAL, COMPLEX } csr_value_type_t;

typedef struct crs {
  /* valide uniquement pour des matrices en double !! */
  csr_value_type_t VALUE_TYPE; 
  int M,N, NZNUM;
  int * COL_IND; 
  int * ROW_PTR;
  double *VAL;
} crs_t;

double crsget(struct crs * spm , int row, int col) ;
double crsgetC(struct crs * spm , int row, int col) ;
void prettyCrs(char * name, struct crs *spm) ;
int CrsIsSymmetric(struct crs * spm);
// multiplication matrice creuse x vecteur dense
//int slowCrsMVM (struct crs * A, double* V, double *RES) ;
//int fastCrsMVM (struct crs * A, double* x, double *y) ;
//
//int vfastCrsMVM (int precision, struct crs * A,
//		 VPFloatArray x, VPFloatArray y);

//double residualNorm2(struct crs * A, double *b, double * x) ;
//int readMMfile(char * filename, struct crs * SPM ) ;
void dumpCrs(struct crs * spm) ;
void crs2dense(struct crs *spm, double * M);
int crsIdent(int n, struct crs * Ident) ;
void crsFree(struct crs * spm);

#ifdef __cplusplus
}
#endif

#endif
