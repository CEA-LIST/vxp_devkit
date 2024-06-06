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
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <time.h>


#include <math.h>

#include "VPSDK/VPFloat.hpp"
#include "VPSDK/VMath.hpp"

#include "data/random_1001.h" 

using namespace VPFloatPackage;

#ifdef __riscv

void * __dso_handle = NULL;

#endif /* __riscv */

/* 
 * les paramètres
 * ==========================
 */
// TOLERANCE PAR DEFAUT: les iterations stoppent quand ||r|| < TOLERANCE
// pratiquement, configurer la tolerance  avec l'option -t 1e-6 par ex
// valeurs raisonables 1e-6<<TOLERANCE<<1e-11
#define TOLERANCE 1e-6
// N_PTS valeur par defaut definit le nombre de points de test
//(de precisions differentes)
// configurer ace -p 2 par ex
#define N_PTS 1

#define MAX_N_PTS 12
int PRECS[MAX_N_PTS]={ 53, 64, 96, 128, 142, 144, 160, 192, 256, 308, 412, 512};
#define EXPONENT_SIZE_BITS 7
#define STRIDE_SIZE_BITS 1


int test_vsqrt(VPF512_T * myVPF, int precision) {
    int nb_errors=0;

    double u=pow(2.0,-precision);

    // faisons simple pour commencer
    VPFloatComputingEnvironment::set_precision(precision); 
    VPFloatComputingEnvironment::set_rounding_mode(VPFloatRoundingMode::VP_RNE); 

    short myBis=precision+EXPONENT_SIZE_BITS+1;
    VPFloat abis (EXPONENT_SIZE_BITS, myBis, STRIDE_SIZE_BITS );
    VPFloat y    (EXPONENT_SIZE_BITS, myBis, STRIDE_SIZE_BITS );
    VPFloat delta(EXPONENT_SIZE_BITS, myBis, STRIDE_SIZE_BITS );

    VPFloat a    (myVPF->chunks[0], myVPF->expo, false, EXPONENT_SIZE_BITS, myBis, STRIDE_SIZE_BITS );

    int nchunks=ceil(((1.0*precision)-1.0)/64.0); 

#ifdef DEBUG
    std::cout<<"nchunks : "<<nchunks<<" \n";
#endif

    for (int jj=1; jj<nchunks; jj++) {
        a.mantissaChunk(jj,myVPF->chunks[jj]);
    }
    /* -------------------------------
        * on a notre valeur a
        * calculer y=sqrt(a)
        * ------------------------------- */
    y = VMath::vsqrt(a);

#ifdef DEBUG
    std::cout<<" reciprocal sqrt of "<<double(a)<<"\n is: "<<double(y)<<" \n";
#endif

    /* -------------------------------
        * delta = y-a^2
        * verifier delta < 2 u 
        * --------------------------------- */
    delta=a-(y*y); if ((double)delta <0.0) delta=-delta;
    int faktor=ceil(double(delta/u));

    if (faktor>2) {
        nb_errors+=1;
        std::cerr<<" large error delta = " <<faktor<<" * u \n";
    }

#ifdef DEBUG
    printf(" delta = a-y^2 = %e < %d \n",(double)delta, faktor);
#endif

    return nb_errors;
}

int do_test_for_precision(int precision) {
        int nb_errors=0;

        printf("precision used: %d\n", precision);

        int test_count = 0;

        // Iterate over predifined values
        for (int ii=0; ii<NBVAL; ii++) {

#ifdef DEBUG
                double valapprox=(1.0+(DATEI[ii].chunks[0]*pow(2.0,-64)))*pow(2.0,DATEI[ii].expo);
                printf("elt # %d : val(approx)=%e\n",ii,valapprox);
#endif

            nb_errors += test_vsqrt(&DATEI[ii], precision);

            test_count++;
        }

        printf("  %d errors on %d tests\n",nb_errors, test_count);

        return nb_errors;
}

/*
 * ---------------- MAIN ---------------------
 * c'est pas trop tot !
 */
int main(int argc, char *argv[])
{
    //clock_t t;
    int nb_errors = 0;
    /*
    */
    int opt;
    int n_pts = N_PTS; // nb de precisions a considerer
    int precision=0;
    char resFile[64]="";
  
    while((opt = getopt(argc, argv, "p:ho:n:")) != -1) { 
        switch(opt) { 
        case 'h': 
            exit(-1);
            break; 
        case 'n':
            sscanf(optarg, "%d", &n_pts);
            if (n_pts > MAX_N_PTS) n_pts=MAX_N_PTS;
            std::cerr<<">> nbre de pts (diff. precisions):: "<<n_pts<<"\n";
            break;
        case 'p':
            sscanf(optarg, "%d", &precision);
            if (n_pts > 1024) n_pts=1024;
            std::cerr<<">> precision:: "<<precision<<"\n";
            break;
        case 'o':
            // output vector file
            strcpy(resFile,optarg);
            break; 
        case '?':
            exit(-1);
            break; 
        } 
    } 

    if ( precision == 0 ) {

        // Iterate over defined precisions
        for (int l_precision_index = 0 ; l_precision_index < MAX_N_PTS; l_precision_index++) {

            nb_errors += do_test_for_precision(PRECS[l_precision_index]);

        }

    } else {

            nb_errors += do_test_for_precision(precision);

    }

    if (nb_errors == 0) {
        std::cout << "SUCCESS" << std::endl;
    } else {
        std::cout << "FAIL " << nb_errors << " errors found." << std::endl;
    }

    return(nb_errors);
}
