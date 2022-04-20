#!/bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
set -x                          # Output commands
set -e                          # Abort on errors



################################################################################
# Search
################################################################################

if [ -z "${GSL_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "GSL selected, but GSL_DIR not set.  Checking some places..."
    echo "END MESSAGE"
    
    FILES="include/gsl/gsl_math.h"
    DIRS="/usr /usr/local /usr/local/gsl /usr/local/packages/gsl /usr/local/apps/gsl ${HOME} c:/packages/gsl"
    for dir in $DIRS; do
        GSL_DIR="$dir"
        for file in $FILES; do
            if [ ! -r "$dir/$file" ]; then
                unset GSL_DIR
                break
            fi
        done
        if [ -n "$GSL_DIR" ]; then
            break
        fi
    done
    
    if [ -z "$GSL_DIR" ]; then
        echo "BEGIN MESSAGE"
        echo "GSL not found"
        echo "END MESSAGE"
    else
        echo "BEGIN MESSAGE"
        echo "Found GSL in ${GSL_DIR}"
        echo "END MESSAGE"
    fi
fi



################################################################################
# Build
################################################################################

if [ -z "${GSL_DIR}" -o "${GSL_DIR}" = 'BUILD' ]; then
    echo "BEGIN MESSAGE"
    echo "Building GSL..."
    echo "END MESSAGE"
    
    # Set locations
    THORN=GSL
    NAME=gsl-1.14
    SRCDIR=$(dirname $0)
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    GSL_DIR=${INSTALL_DIR}
    
(
    exec >&2                    # Redirect stdout to stderr
    set -x                      # Output commands
    set -e                      # Abort on errors
    cd ${SCRATCH_BUILD}
    if [ -e ${DONE_FILE} -a ${DONE_FILE} -nt ${SRCDIR}/dist/${NAME}.tar.gz \
                         -a ${DONE_FILE} -nt ${SRCDIR}/GSL.sh ]
    then
        echo "GSL: The enclosed GSL library has already been built; doing nothing"
    else
        echo "GSL: Building enclosed GSL library"
        
        # Should we use gmake or make?
        MAKE=$(gmake --help > /dev/null 2>&1 && echo gmake || echo make)
        # Should we use gtar or tar?
        TAR=$(gtar --help > /dev/null 2> /dev/null && echo gtar || echo tar)
        
        # Set up environment
        unset LIBS
        if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
            export OBJECT_MODE=64
        fi
        
        echo "GSL: Preparing directory structure..."
        mkdir build external done 2> /dev/null || true
        rm -rf ${BUILD_DIR} ${INSTALL_DIR}
        mkdir ${BUILD_DIR} ${INSTALL_DIR}
        
        echo "GSL: Unpacking archive..."
        pushd ${BUILD_DIR}
        ${TAR} xzf ${SRCDIR}/dist/${NAME}.tar.gz
        
        echo "GSL: Configuring..."
        cd ${NAME}
        ./configure --prefix=${GSL_DIR}
        
        echo "GSL: Building..."
        ${MAKE}
        
        echo "GSL: Installing..."
        ${MAKE} install
        popd
        
        echo "GSL: Cleaning up..."
        rm -rf ${BUILD_DIR}
        
        date > ${DONE_FILE}
        echo "GSL: Done."
    fi
)

    if (( $? )); then
        echo 'BEGIN ERROR'
        echo 'Error while building GSL. Aborting.'
        echo 'END ERROR'
        exit 1
    fi

fi



################################################################################
# Configure Cactus
################################################################################

# Set options
if [ -x ${GSL_DIR}/bin/gsl-config ]; then
    GSL_INC_DIRS=`${GSL_DIR}/bin/gsl-config --cflags | sed -e 's/ \+-[^I][^ ]\+//g;s/^ *-[^I][^ ]\+ *//g;s/ \+-I/ /g;s/^ *-I//g'`;
    GSL_LIB_DIRS=`${GSL_DIR}/bin/gsl-config --libs   | sed -e 's/ \+-[^L][^ ]\+//g;s/^ *-[^L][^ ]\+ *//g;s/ \+-L/ /g;s/^ *-L//g'`;
    GSL_LIBS=`${GSL_DIR}/bin/gsl-config --libs       | sed -e 's/ \+-[^l][^ ]\+//g;s/^ *-[^l][^ ]\+ *//g;s/ \+-l/ /g;s/^ *-l//g'`;
fi

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "HAVE_GSL     = 1"
echo "GSL_DIR      = ${GSL_DIR}"
echo "GSL_INC_DIRS = ${GSL_INC_DIRS}"
echo "GSL_LIB_DIRS = ${GSL_LIB_DIRS}"
echo "GSL_LIBS     = ${GSL_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(GSL_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(GSL_LIB_DIRS)'
echo 'LIBRARY           $(GSL_LIBS)'
