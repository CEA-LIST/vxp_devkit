
mkdir Downloads Extract Install

# GMP Installation
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz -O ./Downloads/gmp-6.3.0.tar.xz
tar -xvf ./Downloads/gmp-6.3.0.tar.xz --directory ./Extract/
cd ./Extract/gmp-6.3.0
./configure --prefix=`readlink -f ../../Install`
make -j 4
make install
cd ../..

# MPFR installation
wget https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.gz -O ./Downloads/mpfr-4.1.0.tar.gz
tar -zxvf ./Downloads/mpfr-4.1.0.tar.gz --directory ./Extract/
cd ./Extract/mpfr-4.1.0
./configure --with-gmp=`readlink -f ../../Install` --prefix=`readlink -f ../../Install` 
make -j 4
make -j 4 check
make install
cd ../..

# MPC installation
wget https://ftp.gnu.org/gnu/mpc/mpc-1.3.1.tar.gz -O ./Downloads/mpc-1.3.1.tar.gz
tar -zxvf ./Downloads/mpc-1.3.1.tar.gz --directory ./Extract/ 
cd ./Extract/mpc-1.3.1
./configure --with-mpfr=`readlink -f ../../Install` --with-gmp=`readlink -f ../../Install` --prefix=`readlink -f ../../Install`
make -j 4
make install
cd ../..

# OSKI installation
wget https://sourceforge.net/projects/oski/files/latest/download -O ./Downloads/oski-1.0.1h.tar.bz2
tar -jxvf  ./Downloads/oski-1.0.1h.tar.bz2 --directory ./Extract/
cd ./Extract/oski-1.0.1h
CFLAGS="-O3 -malign-double -std=c99"  ./configure --enable-int-double --enable-int-dcomplex --prefix=`readlink -f ../../Install`
make -j 4
make -j 4 benchmarks
make install 
cd ../..

rm -Rf Dowloads Extract

mkdir -p Install/pkgconfig

for PC_FILE in ./pkgconfig/*
do
    sed -e "s%@@@INSTALL_PREFIX@@@%`readlink -f ./Install/`%g" $PC_FILE > ./Install/pkgconfig/`basename ${PC_FILE}`
done


