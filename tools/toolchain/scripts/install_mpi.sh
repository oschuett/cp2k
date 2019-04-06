#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${BUILDDIR}"/toolchain.conf

# MPI libraries
case "$MPI_MODE" in
    mpich)
        "${SCRIPTDIR}"/install_mpich.sh "${cp2k_with_mpich}"; load "${BUILDDIR}/setup_mpich"
        ;;
    openmpi)
        "${SCRIPTDIR}"/install_openmpi.sh "${cp2k_with_openmpi}"; load "${BUILDDIR}/setup_openmpi"
        ;;
esac

# update toolchain config again
export -p > ${BUILDDIR}/toolchain.conf
