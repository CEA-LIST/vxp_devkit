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
 *  @file        vblas_complex.c
 *  @author      Jerome Fereyre
 */

#include "VRPSDK/vblas.h"

void cscal (    int n, void * alpha , void *x, vpfloat_off_t x_inc ) {

}

void csscal (   int n, float alpha , void *x, vpfloat_off_t x_inc ) {

}

void zscal (    int n, void * alpha , void *x, vpfloat_off_t x_inc ) {

}

void zsscal (   int n, double alpha , void *x, vpfloat_off_t x_inc ) {

}
void vcscal (   int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                void *x, vpfloat_evp_t x_evp ) {

}
void vcsscal (  int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                void *x, vpfloat_evp_t x_evp ) {
                
}  

void ccopy (    int n,
                const void *x, vpfloat_off_t x_inc ,
                void *y, vpfloat_off_t y_inc ) {

}

void zcopy (    int n,
                const void *x, vpfloat_off_t x_inc ,
                void *y, vpfloat_off_t y_inc ) {

}

void vccopy (   int n,
                const void *x, vpfloat_evp_t x_evp ,
                void *y, vpfloat_evp_t y_evp ) {

}

void caxpy (int n, void * alpha ,
            const void *x, vpfloat_off_t x_inc ,
            void *y, vpfloat_off_t y_inc ) {

}

void zaxpy (int n, void * alpha ,
            const void *x, vpfloat_off_t x_inc ,
            void *y, vpfloat_off_t y_inc ) {

}

void vcaxpy (   int precision , int n,
                const void *alpha , vpfloat_evp_t a_evp ,
                const void *x, vpfloat_evp_t x_evp ,
                void *y, vpfloat_evp_t y_evp ) {

}

void cdotu (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r) {
}

void cdotc (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r) {
}

void zdotu (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r) {
}

void zdotc (int n,
            const void *x, vpfloat_off_t x_inc,
            const void *y, vpfloat_off_t y_inc,
            void * r) {
}

void vcdotu (   int precision , int n,
                const void *x, vpfloat_evp_t x_evp ,
                const void *y, vpfloat_evp_t y_evp ,
                void *r, vpfloat_evp_t r_evp ) {

}

void vcdotc (   int precision , int n,
                const void *x, vpfloat_evp_t x_evp ,
                const void *y, vpfloat_evp_t y_evp ,
                void *r, vpfloat_evp_t r_evp ) {

}                

void cgemv (char trans , int m, int n,
            const void * alpha ,
            const void *a, int lda ,
            const void *x, int x_inc ,
            const void * beta ,
            void *y, int y_inc ) {

}

void zgemv (char trans , int m, int n,
            const void * alpha ,
            const void *a, int lda ,
            const void *x, int x_inc ,
            const void * beta ,
            void *y, int y_inc ) {

}

void vcgemv (   int precision , char trans , int m, int n,
                const void *alpha , vpfloat_evp_t alpha_evp ,
                const void *a, vpfloat_evp_t a_evp , int lda, 
                const void *x, vpfloat_evp_t x_evp ,
                const void *beta , vpfloat_evp_t beta_evp ,
                void *y, vpfloat_evp_t y_evp ) {

}