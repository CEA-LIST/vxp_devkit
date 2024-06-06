# Copyright 2022 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##
#  author:      Cesar Fuguet-Tortolero
#  date:        November 18, 2022
#  description: Script to cross-compile the GMP and MPFR libraries
##
#!/bin/bash

THIS_SCRIPT=$0
THIS_DIRECTORY=$(dirname $(readlink -f ${THIS_SCRIPT}))
BUILD_GMP=1
BUILD_MPFR=1
BUILD_RISCV=1
BUILD_LINUX=1

if [ -z ${INSTALL_ROOT} ]
then
	INSTALL_ROOT="${PWD}"
fi

if [ -z ${LINUX_VERSION} ]
then
	LINUX_VERSION="-local"
fi


# set common environment settings
CC=${RISCV_HOME}/bin/riscv64-unknown-elf-gcc
LOG_DIRECTORY=${THIS_DIRECTORY}/build

build_gmp()
{
    CC=${1}
    CFLAGS=${2}
    LDFLAGS=${3}
    BUILD_DIRECTORY=${4}
    INSTALL_DIRECTORY=${5}
    SRC_DIRECTORY=${6}
    HOST_OPTION=${7}

    (
    mkdir -p ${INSTALL_DIRECTORY}

    mkdir -p ${BUILD_DIRECTORY}/gmp
    cd ${BUILD_DIRECTORY}/gmp

    # configure
    ${SRC_DIRECTORY}/gmp-6.2.0/configure \
        CC=${CC} \
        CFLAGS="${CFLAGS}" \
        LDFLAGS="${LDFLAGS}" \
        ${HOST_OPTION} \
        --prefix=${INSTALL_DIRECTORY}

    # build
    make -j3
    if [[ $? != 0 ]] ; then
        echo "error: build of the gmp library failed"
        exit 1 ;
    fi

    # install
    make install
    )
}


build_mpfr()
{
    CC=${1}
    CFLAGS=${2}
    LDFLAGS=${3}
    BUILD_DIRECTORY=${4}
    INSTALL_DIRECTORY=${5}
    SRC_DIRECTORY=${6}
    HOST_OPTION=${7}

    (
    mkdir -p ${INSTALL_DIRECTORY}

    mkdir -p ${BUILD_DIRECTORY}/mpfr
    cd ${BUILD_DIRECTORY}/mpfr

    # configure
    ${SRC_DIRECTORY}/mpfr-4.1.0/configure \
        CC=${CC} \
        CFLAGS="${CFLAGS}" \
        LDFLAGS="${LDFLAGS}" \
        --with-gmp=${INSTALL_DIRECTORY} \
        ${HOST_OPTION} \
        --prefix=${INSTALL_DIRECTORY}

    # build
    make -j3
    if [[ $? != 0 ]] ; then
        echo "error: build of the gmp library failed"
        exit 1 ;
    fi

    # install
    make install ;
    )
}

mkdir -p archives
if [[ ! -e archives/mpfr-4.1.0.tar.xz ]] ; then
    curl -o archives/mpfr-4.1.0.tar.xz https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.xz
fi
if [[ ! -e archives/gmp-6.2.0.tar.xz ]] ; then
    curl -o archives/gmp-6.2.0.tar.xz https://ftp.gnu.org/gnu/gmp/gmp-6.2.0.tar.xz
fi

mkdir -p src/
[[ ! -d src/gmp-6.2.0 ]] && tar xvf archives/gmp-6.2.0.tar.xz -C src/ ;
[[ ! -d src/mpfr-4.1.0 ]] && tar xvf archives/mpfr-4.1.0.tar.xz -C src/ ;

if [ ${BUILD_RISCV} ]
then

	CFLAGS='-mabi=lp64d -march=rv64imafdc -mcmodel=medany'
	LDFLAGS='-mabi=lp64d -march=rv64imafdc -mcmodel=medany'
	HOST_OPTION='--host=riscv64'
	BUILD_DIRECTORY="${PWD}/build-hard-float"
	INSTALL_DIRECTORY="${INSTALL_ROOT}/riscv-hard-float"

	# build the GMP library
	if [[ ${BUILD_GMP} == 1 ]] ; then
		build_gmp "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi

	# build the MPFR library
	if [[ ${BUILD_MPFR} == 1 ]] ; then
		build_mpfr "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi

	CFLAGS='-mabi=lp64 -march=rv64imac -mcmodel=medany'
	LDFLAGS='-mabi=lp64 -march=rv64imac -mcmodel=medany'
	HOST_OPTION='--host=riscv64'
	BUILD_DIRECTORY="${PWD}/build-soft-float"
	INSTALL_DIRECTORY="${INSTALL_ROOT}/riscv-soft-float"

		# build the GMP library
	if [[ ${BUILD_GMP} == 1 ]] ; then
		build_gmp "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi

	# build the MPFR library
	if [[ ${BUILD_MPFR} == 1 ]] ; then
		build_mpfr "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi

fi


if [ ${BUILD_LINUX} ]
then
	CC=/usr/bin/gcc
	CFLAGS=''
	LDFLAGS=''
	HOST_OPTION=''
	BUILD_DIRECTORY="${PWD}/build-linux"
	INSTALL_DIRECTORY="${INSTALL_ROOT}/x86_64-linux${LINUX_VERSION}"

	# build the GMP library
	if [[ ${BUILD_GMP} == 1 ]] ; then
		build_gmp "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi

	# build the MPFR library
	if [[ ${BUILD_MPFR} == 1 ]] ; then
		build_mpfr "${CC}" "${CFLAGS}" "${LDFLAGS}" "${BUILD_DIRECTORY}" "${INSTALL_DIRECTORY}" "${PWD}/src" "${HOST_OPTION}"
	fi
fi
