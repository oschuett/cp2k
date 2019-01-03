#!/bin/bash -e

# author: Ole Schuett

if (( $# < 2 )) || (( 3 < $# )) ; then
    echo "Usage: test_regtest.sh <ARCH> <VERSION> [ <TESTOPTS> ]"
    exit 1
fi

ARCH=$1
VERSION=$2
TESTOPTS="${TESTOPTS} $3"

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Running Regtests =========="
cd /workspace/cp2k
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

#EOF
