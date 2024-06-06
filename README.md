
# How to use the VRP development kit

## Clone the git repository containing the VRP development kit

```
git clone  --branch v0.1-TRISTAN  ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/epi/vrp-software-dev-env.git v0.1-TRISTAN

cd v0.1-TRISTAN

git submodule update --init --recursive
```
## Compilation of external libraries used by the VRP development kit

This step may take a long time depending on the performance of your system (compilation of several software libraries)

```
cd external_libraries

bash ./install.sh
```

## Compiling the VRP development kit

To compile the VRP development kit you have to run the build.sh script with the following options : --rebuild --all
This version of the development kit supports only MPFR software emulation of variable precision.

```
bash ./build.sh --rebuild --all
```

### Using the VRP development kit for a new program

To use the VRP development kit you need to update the content of several environment variables : PATH, LD_LIBRARY_PATH, PKG_CONFIG_PATH.

To simplify these modification it is possible to use the "module" tool for which we provide a set of configuration files in the "modules" directory

To do so you have to:

1. Indicate to the module tool where the configuration files are located (the "module" command must be run from the cloned directory root)

```
module use ./modules
```

2. Load the main module that sets the appropriate environment

The "module" tool manage dependencies, so you have to load only one module by hand, the other mandatory modules will be loaded automatically.

```
module load vrp/sdk/linux_x86_64/current-dev-linux_x86_64
```


### Launch VRP development kit non-regression tests

There are some non-regression tests allowing to check that VRP development kit works properly.
You can find them in the vp_sdk/vp_sdk_unit_tests/ directory.
To run them you have to setup environment variables with the module command (previous step) and run the nonreg.sh script

```
./nonreg.sh --linux --exclude eigen --exclude vsqrt --exclude matrix
```

# Create a new application using the VRP development kit

1. You shall build the environment as described in previous sections.

2. You shall setup environment variables (using the "module" command).

3. Create a Makefile to build your application. Use the following example as skeleton, it contains the pkg-config command example needed to compile and link your application with on the VRP development kit.

```
TARGETS=test
BUILD_DIR=$(shell readlink -f ./build)

CXXFLAGS=$(shell pkg-config --cflags vp_sdk_linux_x86_64 matrix_sdk_linux_x86_64) -ggdb -O0 -Wall
LDFLAGS=$(shell pkg-config --libs vp_sdk_linux_x86_64 matrix_sdk_linux_x86_64)

all: clean ${TARGETS}

clean:
	-rm -Rf $(BUILD_DIR) test

test: ${BUILD_DIR}/test.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) -lm

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
```

4. Code your application. Here is an example. The Makefile in the previous step considers that the following application is in a test.c file. You can use another name but remember to modify the Makefile accordingly.

```
#include <iostream>
#include <OSKIHelper.hpp>
#include <Matrix/matrix.h>
#include "VPSDK/VBLAS.hpp"
#include "VRPSDK/perfcounters/cpu.h"

using namespace VPFloatPackage;

#define ITER_MAX 5

void initB(double * B, int n) {
    for ( int index = 0 ; index < n ; index++ ) {
        B[index] = 1.0;
    }
}

int my_cg(  int precision,
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

    /* p_k <- b */
    VBLAS::vcopy(n,b,p_k);
    /* x_k = {0} */
    VBLAS::vzero(precision, n, x_k);
    VBLAS::vzero(precision, n, Ap_k);
    rs_next = 0.0;

    // rs = rk'*rk
    VBLAS::vdot(precision, n, r_k, r_k, rs);

    for (nbiter = 0; nbiter < n*ITER_MAX; ++nbiter) {
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
    
        // alpha = Ap_k * p_k
        VBLAS::vdot(precision, n, p_k, Ap_k,  alpha);

        // alpha = rs/alpha
        alpha = rs/alpha;

        // x_k = x_k + alpha*p_k
        VBLAS::vaxpy(precision, n, alpha, p_k,  x_k);
	
        // r_k = r_k - alpha * Ap_k
        VBLAS::vaxpy(precision, n, -alpha, Ap_k, r_k);


        // rs_next (rs:r square)
        VBLAS::vdot(precision, n, r_k,  r_k, rs_next);

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

int main(int argc, char ** argv) {

    /* Loading a matrix from matrix market repo and convert it into matrix_t type */
    char * l_matrix_path = argv[1];

    std::cout << "Loading input matrix: " << l_matrix_path << "..." << std::endl;

    oski_matrix_wrapper_t l_my_oski_csr_matrix = OSKIHelper::loadFromFile(l_matrix_path);

    matrix_t l_my_csr_matrix = OSKIHelper::toMatrix(l_my_oski_csr_matrix);

    matrix_t l_my_dense_matrix = OSKIHelper::toDense(l_my_oski_csr_matrix);

    std::cout << "matrix loaded!" << std::endl;

    double * X = (double *)malloc(sizeof(double) * l_my_csr_matrix->n);
    memset(X, 0, sizeof(double) * l_my_csr_matrix->n);
    double * B = (double *)malloc(sizeof(double) * l_my_csr_matrix->n);
    memset(B, 0, sizeof(double) * l_my_csr_matrix->n);

    initB(B, l_my_csr_matrix->n);

    VPFloatArray Xv(X, l_my_csr_matrix->n);
    VPFloatArray Bv(B, l_my_csr_matrix->n);

    int l_iteration_count_dense = my_cg(    500,
                                            0, // Trans
                                            l_my_csr_matrix->n, //n
                                            Xv,  // valeur de sortie 
                                            l_my_dense_matrix,
                                            Bv,
                                            1e-12,
                                            11,
                                            1);

    int l_iteration_count_sparse = my_cg(   500,
                                            0, // Trans
                                            l_my_csr_matrix->n, //n
                                            Xv,  // valeur de sortie 
                                            l_my_csr_matrix,
                                            Bv,
                                            1e-12,
                                            11,
                                            1);

    std::cout << "Iteration : sparse=" << l_iteration_count_sparse <<  " - dense=" << l_iteration_count_dense << std::endl;

    return 0;
}
```

5. Build your application with the "make" command

```
make clean
make
```

6. Run your application. For the previous code example you have to get a matrix in "mtx" file format. You can find an example on nist web site. Download an example :

```
wget https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc3/bcsstk22.mtx.gz && gunzip bcsstk22.mtx.gz
```

Run the test program on this matrix:

```
./test bcsstk22.mtx

Loading input matrix: bcsstk22.mtx...
 symmetric matrix 138 x 138 with 417 non-zeros
Set VALUE_TYPE to REAL
LDA set automatically to 152
matrix loaded!
Iteration : sparse=157 - dense=157
```
