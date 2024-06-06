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
#include "MTXUtil/mmio.h"

int printVectorFile(char * filename, int N, double * VEKTOR ) {
  FILE *f;
  int ii;
  if ((f = fopen(filename, "w")) == NULL) {
    fprintf(stderr,"could not write vector file %s\n",filename);
    exit(-1);
  } else {
    for (ii=0; ii<N; ii++) {
      fprintf(f,"%e \n",VEKTOR[ii]);
    }
  }
  return(N);
 }
int printFFVectorFile(char * filename, int N, double * VEKTOR ) {
  // specific Freefem++
  FILE *f;
  int ii;
  if ((f = fopen(filename, "w")) == NULL) {
    fprintf(stderr,"could not write vector file %s\n",filename);
    exit(-1);
  } else {
    fprintf(f,"%d \n",N);
    for (ii=0; ii<N; ii++) {
      fprintf(f,"%e \n",VEKTOR[ii]);
    }
  }
  return(N);
 }

int readVectorFile(char * filename, int N, double * VEKTOR ) {
  // lit un fichier 
  FILE *f;
  int ii;

  /*
   *  1: open & parse file
   */
  fprintf(stderr,"open & parse file %s\n",filename);

  if ((f = fopen(filename, "r")) == NULL) {
    fprintf(stderr,"could not find vector file %s\n",filename);
    exit(-1);
  } else {
    for (ii=0; ii<N; ii++) {
      if (fscanf(f, "%le\n", &VEKTOR[ii]) == 0) exit(-1);
      //fprintf(stderr,"je vois VEKTOR[%ii]=%e\n",ii,VEKTOR[ii]);
    }
  }
  return(N);
}

int readMMfile(char * filename, crs_t * SPM ) {
  // lit le fchier matrix market
  // renseigne la structure crs spm
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;   
  int i, *ROW, *COL;
  double *val;
  /*
   *  1: open & parse file
   */
  fprintf(stderr,"open & parse file %s\n",filename);

  if ((f = fopen(filename, "r")) == NULL) {
    fprintf(stderr,"could not find this file %s\n",filename);
    exit(-1);
  }

  if (mm_read_banner(f, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }


  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
      mm_is_sparse(matcode) ) {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
      exit(1);
  fprintf(stderr,"// ");
  if (mm_is_symmetric(matcode)) fprintf(stderr, " symmetric ");
  fprintf(stderr,"matrix %d x %d with %d non-zeros\n",M,N,nz);

  /* reserve memory for matrices */

  if (mm_is_symmetric(matcode)) {
    ROW = (int *) malloc(2* nz * sizeof(int));
    COL = (int *) malloc(2* nz * sizeof(int));
    val = (double *) malloc(2*nz * sizeof(double));
  } else {
    ROW = (int *) malloc(nz * sizeof(int));
    COL = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
  }

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
  int nzidx = nz;

  for (i=0; i<nz; i++)
  {
    if (fscanf(f, "%d %d %lg\n", &ROW[i], &COL[i], &val[i]) == 0) exit(-1);
	  //fprintf(stderr,"j ai vu ligne %d:: %d %d %f\n",i,ROW[i],COL[i],val[i]);

    if (mm_is_symmetric(matcode)) {
      // je rajoute le symetrique a la fin
      if (ROW[i]!=COL[i]) {
        ROW[nzidx]=COL[i];
        COL[nzidx]=ROW[i];
        val[nzidx]=val[i];
        nzidx++;
      }
    }
  }
  nz=nzidx; 

  if (f !=stdin) fclose(f);

  // debug
  //for (int nn=0; nn < nzidx; nn++)
  //  printf("%d:: %d %d %e \n",nn,ROW[nn],COL[nn],val[nn]);
  /*
    * 2: order ROW & COL 
    */
  int *ROW_ORD, *COL_ORD;
  double *val_ORD;
  ROW_ORD = (int *) malloc(nz * sizeof(int));
  COL_ORD = (int *) malloc(nz * sizeof(int));
  val_ORD = (double *) malloc(nz * sizeof(double));
  int workingIdx, col;
  int mm, kk, pROW;
  int *curRowIdx;
  curRowIdx = (int *) malloc((nz+1) * sizeof(int));
  workingIdx=0; // pointeur de travail pour remplir ROM_ORD/COL_ORD/val_ORD
  for (mm=1; mm<=M; mm++) {
    // invalider les cases a priori
    for (kk=0; kk<=nz; kk++) curRowIdx[kk]=-1; 
      // prendre les lignes 1 par 1 avec mm (mm ds [1,M])
      // parcours de ROW
      for (pROW=0; pROW<nz; pROW++) {
        if (ROW[pROW]==mm) {
          // col=COL[pROW] pointe sur la colonne ds [1:N]
          col=COL[pROW]; // dans [1:N]
          // on stocke pRoW dans curRowIdx[col]
          curRowIdx[col]=pROW;
        }
    }
    // dans curRowIdx on a tout le rang mm, reste a ranger dans ROW_ORD, etc.
    //printf("\n curRowIdx[%d]: \n",mm);
    //for(int tt=1; tt<nz+1; tt++) printf(" %d ",curRowIdx[tt]);
    
    for (kk=1; kk<=nz; kk++) {
      if (curRowIdx[kk]>=0) { // sinon c'est vide
        // curRowIdx[kk] pointe sur l'element suivant, qu'on va ranger
        // dans ROW_ORD[workingIdx]
        ROW_ORD[workingIdx]=mm;
        COL_ORD[workingIdx]=COL[curRowIdx[kk]];
        val_ORD[workingIdx]=val[curRowIdx[kk]];
        workingIdx++;
      }; 
    }
  }

  free(curRowIdx);

  // a la fin on devrait avoir idx=nz
  if (workingIdx!=nz) fprintf(stderr,"mauvais decompte pendant la phase de tri\n");
  // pas oublier de liberer ROW & COL & val a la fin
  free(ROW);
  free(COL);
  free(val);

  /* 
    ***********************
    * 3: create my crs struct 
    *
    */

  //printf("initialisation\n");
  SPM->M=M;
  SPM->N=N;
  SPM->NZNUM=nz;
  SPM->VAL=val_ORD;
  // construire COL_IND & ROW_PTR
  // les colonnes indexees de 0 à _M-1
  SPM->COL_IND =(int *) malloc(SPM->NZNUM*sizeof(int));
  SPM->ROW_PTR =(int *) malloc((SPM->M+1)*sizeof(int));
  
  int idx;       // parcourt la triple liste _VAL, _ROW, _COL
  int new_row=1, previous_row=0; // logic
  SPM->ROW_PTR[0]=ROW_ORD[0];

  for (idx=0; idx<SPM->NZNUM; idx++) {
    SPM->COL_IND[idx]=COL_ORD[idx];
    // printf("je vois _ROW[%d]=%d",idx,ROW[idx]);
    if (ROW_ORD[idx]==previous_row) {
      // printf(" ... je saute \n");
      // skip
    } else {
      new_row=ROW_ORD[idx];
      // printf(" ... je note ROW_PTR[%d]=%d\n",new_row-1,idx+1);
      // attention ! les valeurs de ROW_PTR sont de 1 à N -> add 1
      // le stockage ds ROW_PTR sont de 0 a M-1 -> retrancher 1 !!

      for (kk=previous_row; kk <= (new_row-1); kk++)
        //SPM->ROW_PTR[new_row-1]=idx+1;
        SPM->ROW_PTR[kk]=idx+1;
      }
      previous_row=new_row;

      SPM->ROW_PTR[M]=nz+1; // convention IMPORTANTE
  }
  
  free(COL_ORD);
  free(ROW_ORD);

  return(ret_code);
}

/*
*/
