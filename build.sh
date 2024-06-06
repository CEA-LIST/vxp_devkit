#!/bin/bash

######################################################################################
#  Setup various PATH needed to retrieve or install libraries

EXTERNAL_LIBS_INSTALL_PATH=`readlink -f external_libraries/Install`
MPFR_LIBRARY_ROOT_PATH=${EXTERNAL_LIBS_INSTALL_PATH}
OSKI_LIBRARY_ROOT_PATH=${EXTERNAL_LIBS_INSTALL_PATH}

MODULE_LIST="matrix_sdk vrp_sdk vp_sdk"
LINUX_AVAILABLE_PLATFORM_LIST="linux_x86_64 "
AVAILABLE_PLATFORM_LIST="${LINUX_AVAILABLE_PLATFORM_LIST}"

declare -A MODULE_LIST_PER_PLATFORM=(
['linux_x86_64']="matrix_sdk vrp_sdk vp_sdk "
)

LINUX_APPS="vrp_run_bare_app vrp_test_solver "


SELECTED_PLATFORM_LIST=""


function do_on_each_module() {
    CMD=${1}
    MODULE_LIST=${2}
    
    for MODULE in ${MODULE_LIST}
    do
        mkdir -p ${MODULE}

        cd ${MODULE}

        eval ${CMD}

        if [ ${?} -ne 0 ]
        then
            exit 1
        fi

        cd ..
    done
}


CLEAN_BEFORE_BUILD=0
BUILD_FOR_LINUX=0
BUILD_TYPE="Debug"

######################################################################################
# Processing command line arguments 
if [ ${#} -ne 0 ] 
then

	while [ ${#} -ne 0 ]
	do
		case ${1} in
			--help)
				echo "--all                     : Build everything"
				echo "--linux                   : Build for Linux platform"
                echo "--rebuild                 : Force rebuild"
				exit 0
				;;
			--all)
				BUILD_FOR_LINUX=1
                SELECTED_PLATFORM_LIST="${AVAILABLE_PLATFORM_LIST}"
				;;
			--linux)
				BUILD_FOR_LINUX=1
                SELECTED_PLATFORM_LIST="${SELECTED_PLATFORM_LIST} linux_x86_64"
				;;
            --rebuild)
                CLEAN_BEFORE_BUILD=1
                ;;
            --release)
                BUILD_TYPE="Release"
                ;;
			*)
				echo "unknown option ${1}"
				exit1
				;;
		esac
		shift
	done
	
else
	BUILD_FOR_LINUX=1
fi


CMAKE_COMMON_OPTIONS=" -j10 -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON "

ROOT_PATH=${PWD}
INSTALL_ROOT_PATH=${ROOT_PATH}/new_install/${BUILD_TYPE}
BUILD_ROOT_PATH=${ROOT_PATH}/new_build/${BUILD_TYPE}
VRP_LINUX_DRIVER_INSTALLATION_PATH="${INSTALL_ROOT_PATH}/linux_x86_64/vrp_linux_driver/"


if [ ${BUILD_FOR_LINUX} -eq 1 ]
then
    mkdir -p ${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include
    if [ ! -f ${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include/vrp_ioctl.h ] || [ external_includes/vrp_ioctl.h -nt ${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include/vrp_ioctl.h ]
    then
        cp external_includes/vrp_ioctl.h ${VRP_LINUX_DRIVER_INSTALLATION_PATH}/include
    fi
fi

mkdir -p ${BUILD_ROOT_PATH} 

cd ${BUILD_ROOT_PATH} 

ORIGINAL_PKG_CONFIG_PATH="${PKG_CONFIG_PATH}"

for PLATFORM in ${SELECTED_PLATFORM_LIST}
do
    # Build the PKG_CONFIG_PATH variable for next platform compilation
    PKG_CONFIG_PATH="${ORIGINAL_PKG_CONFIG_PATH}:${MPFR_LIBRARY_ROOT_PATH}/pkgconfig:${OSKI_LIBRARY_ROOT_PATH}/pkgconfig"

    for MODULE in ${MODULE_LIST_PER_PLATFORM[${PLATFORM}]}
    do
        PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${INSTALL_ROOT_PATH}/${PLATFORM}/${MODULE}/pkgconfig"
    done

    if [ ${CLEAN_BEFORE_BUILD} -eq 1 ]
    then
        rm -Rf ${PLATFORM}
    fi

    # If PLATFORM directory do no exists we need to initialize the build directory
    if [[ ! -d ${PLATFORM} ]]
    then
        mkdir -p ${PLATFORM}

        cd ${PLATFORM}

        PLATFORM_CMAKE_OPTIONS="${CMAKE_COMMON_OPTIONS} -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT_PATH}/${PLATFORM}/\${MODULE}"


        PLATFORM_CMAKE_OPTIONS="${PLATFORM_CMAKE_OPTIONS} -C ${ROOT_PATH}/\${MODULE}/linux_configuration_parameters.cmake"
        PLATFORM_CMAKE_OPTIONS="${PLATFORM_CMAKE_OPTIONS} -DVRP_LINUX_DRIVER_INSTALLATION_PATH=${VRP_LINUX_DRIVER_INSTALLATION_PATH}"

        ##############################################################################
        #
        do_on_each_module "PKG_CONFIG_PATH=${PKG_CONFIG_PATH} cmake ${PLATFORM_CMAKE_OPTIONS} ${ROOT_PATH}/\${MODULE} && cmake --build . && cmake --install ." "${MODULE_LIST_PER_PLATFORM[${PLATFORM}]}"

        cd ..
    else
        cd ${PLATFORM}

        ##################################################################################
        #
        do_on_each_module "cmake --build . && cmake --install ." "${MODULE_LIST_PER_PLATFORM[${PLATFORM}]}"

        cd ..
    fi

done

cd ${ROOT_PATH}

######################################################################################
# Build the linux apps sample application 
if [ ${BUILD_FOR_LINUX} -eq 1 ]
then
    # Build the PKG_CONFIG_PATH variable for next platform compilation
    PKG_CONFIG_PATH="${ORIGINAL_PKG_CONFIG_PATH}:${MPFR_LIBRARY_ROOT_PATH}/pkgconfig:${OSKI_LIBRARY_ROOT_PATH}/pkgconfig"

    for MODULE in ${MODULE_LIST_PER_PLATFORM['linux_x86_64']}
    do
        PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${INSTALL_ROOT_PATH}/${PLATFORM}/${MODULE}/pkgconfig"
    done

	cd vrp_sdk_linux_apps

    for CURRENT_APP in ${LINUX_APPS}
    do
        cd ${CURRENT_APP}

        export PKG_CONFIG_PATH

        if [ ${CLEAN_BEFORE_BUILD} -eq 1 ]
        then
            make -f makefile.vrp_sdk clean 
        fi

        make -f makefile.vrp_sdk 

        if [ ${?} -ne 0 ]
        then
            exit 1
        fi

        mkdir -p ${INSTALL_ROOT_PATH}/linux_x86_64/bin

        cp ${CURRENT_APP} ${INSTALL_ROOT_PATH}/linux_x86_64/bin

        cd ..
    done

	cd ..
fi
